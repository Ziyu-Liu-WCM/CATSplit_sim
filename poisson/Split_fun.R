singleSplit<-function(rep, inclusion_table, otu_table, taxonomy_table, meta_table, tree_data, 
                      metric,  qval_bound = 0.05){
  
  data1_ind <- sample(colnames(otu_table), ncol(otu_table)/2, replace = FALSE)
  data2_ind <- setdiff(colnames(otu_table), data1_ind)
  
  data1_otu_table <- otu_table[,data1_ind]
  data2_otu_table <- otu_table[,data2_ind]
  
  data1_meta_table <- meta_table[data1_ind,,drop=FALSE]
  data2_meta_table <- meta_table[data2_ind,,drop=FALSE]
  
  testList<-unique(taxonomy_table$taxa_to_genus)
  
  temp1<- computeR2(testList = testList, otutable = data1_otu_table, taxonomy = taxonomy_table, metaData = data1_meta_table,
                     tree = tree_data, metric = metric)
  R21<-temp1[1,1]
  beta1<-temp1[,3]
  beta2 <- computeR2(testList = testList, otutable = data2_otu_table, taxonomy = taxonomy_table, metaData = data2_meta_table,
                     tree = tree_data, metric = metric)[,3]
  R22<-temp1[1,1]
  M <- sign(beta1 * beta2) * (abs(beta1) * abs(beta2))
  selected_index <- analys(M, abs(M), qval_bound)
  M_selected <- M[selected_index]
  inclusion_rep <- ifelse(inclusion_table$feature %in% names(M_selected), 1, 0)
  save(R21,R22,beta1,beta2, file=paste0("savedData/Internal", rep, ".RData"))
  return(inclusion_rep)
}



CATSplit_parallel <- function(otu_table, taxonomy_table, meta_table, tree_data, 
                              metric, nCore, nReps, qval_bound = 0.05){
  
  inclusion_table <- data.frame(feature = unique(taxonomy_table$taxa_to_genus))

  # Register parallel backend
  cl <- makeCluster(nCore)
  registerDoParallel(cl)

  # Data Splitting
  start_time <- Sys.time()

  # Parallel
  results <- foreach(rep = 1:nReps, .combine = 'cbind',
                     .packages = c("ape", "vegan", "GUniFrac", "doParallel"),
                     .export = c("computeR2","permuteRows", "compDist", "analys", "find_tau","singleSplit")) %dopar% {
                       
                       singleSplit(rep,inclusion_table,otu_table, taxonomy_table, meta_table, tree_data, metric, qval_bound)
                     }

  end_time <- Sys.time()
  cat("Run time:", end_time - start_time)
  # Stop the parallel backend
  stopCluster(cl)
  
  # Deregister the parallel backend
  registerDoSEQ()
  
  # Remove objects
  rm(cl)
  
  # Call garbage collection
  gc()
  
  
  ## Calculated the inclusion rate
  column_sums <- ifelse(colSums(results)==0, 1, colSums(results))
  inclusion_table$inclusion_rate <- rowMeans(sweep(results, 2, column_sums, "/"))
  
  inclusion_rate <- inclusion_table$inclusion_rate
  names(inclusion_rate) <- inclusion_table$feature
  inclusion_rate <- inclusion_rate[order(inclusion_rate, decreasing = FALSE)]
  
  
  ## find the largest l
  cum_sum <- cumsum(inclusion_rate)
  l <- max(which(cum_sum <=qval_bound))
  
  ## select features and order them in descending order
  selected_features <- rev(names(inclusion_rate[-(1:l)]))
  
  return(selected_features)
}


CATSplit_noparallel <- function(otu_table, taxonomy_table, meta_table, tree_data, 
                                metric, nReps, qval_bound = 0.05){
  inclusion_table <- data.frame(feature = unique(taxonomy_table$taxa_to_genus))
  
  # Non-parallel
  results<-NULL
  for(iter in 1:nReps){
    results<-cbind(results,singleSplit(iter,inclusion_table,otu_table, taxonomy_table, meta_table, tree_data, metric, qval_bound))
  }

  ## Calculated the inclusion rate
  column_sums <- ifelse(colSums(results)==0, 1, colSums(results))
  inclusion_table$inclusion_rate <- rowMeans(sweep(results, 2, column_sums, "/"))
  
  # write.csv(inclusion_table, file = "C:/Users/lclce/Desktop/2nd Sem/Yushu RA/CATSplit/inclusion_table.csv", row.names = FALSE)
  # inclusion_table <- read.csv("C:/Users/lclce/Desktop/2nd Sem/Yushu RA/CATSplit/inclusion_table.csv")
  
  inclusion_rate <- inclusion_table$inclusion_rate
  names(inclusion_rate) <- inclusion_table$feature
  inclusion_rate <- inclusion_rate[order(inclusion_rate, decreasing = FALSE)]
  
  
  ## find the largest l
  cum_sum <- cumsum(inclusion_rate)
  l <- max(which(cum_sum <= qval_bound))
  
  ## select features and order them in descending order
  selected_features <- rev(names(inclusion_rate[-(1:l)]))
  
  return(selected_features)
}
