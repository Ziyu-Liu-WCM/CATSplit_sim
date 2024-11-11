# metric <- "Weighted UniFrac"
# parallel <- TRUE
# nCore <- 6
# qval_bound <- 0.05
# plot(tree_data)
# summary(tree_data)




CATSplit_parallel <- function(otu_table, taxonomy_table, meta_table, tree_data, 
                              metric, parallel, nCore, nReps, qval_bound = 0.05){
  
  inclusion_table <- data.frame(feature = unique(taxonomy_table$taxa_to_genus))

  # Register parallel backend
  cl <- makeCluster(nCore)
  registerDoParallel(cl)

  # Data Splitting
  start_time <- Sys.time()

  # Parallel
  results <- foreach(rep = 1:nReps, .combine = 'cbind',
                     .packages = c("ape", "vegan", "GUniFrac", "doParallel"),
                     .export = c("computeR2", "compDist", "analys", "find_tau")) %dopar% {
                       ## Data Splitting
                       data1_ind <- sample(colnames(otu_table), ncol(otu_table)/2, replace = FALSE)
                       data2_ind <- setdiff(colnames(otu_table), data1_ind)

                       data1_otu_table <- otu_table[,data1_ind]
                       data2_otu_table <- otu_table[,data2_ind]

                       data1_meta_table <- meta_table[data1_ind,, drop = FALSE]
                       data2_meta_table <- meta_table[data2_ind,, drop = FALSE]

                       testList<-unique(taxonomy_table$taxa_to_genus)

                       beta1 <- computeR2(testList = testList, otutable = data1_otu_table, taxonomy = taxonomy_table, metaData = data1_meta_table,
                                          tree = tree_data, metric = metric, parallel = parallel, nCore = nCore)[,3]
                       beta2 <- computeR2(testList = testList, otutable = data2_otu_table, taxonomy = taxonomy_table, metaData = data2_meta_table,
                                          tree = tree_data, metric = metric, parallel = parallel, nCore = nCore)[,3]

                       M <- sign(beta1 * beta2) * (abs(beta1) * abs(beta2))

                       # tau <- find_tau(M, target_fdr = 0.05)
                       # M_selected <- M[M > tau]
                       
                       selected_index <- analys(M, abs(M), qval_bound)
                       M_selected <- M[selected_index]

                       inclusion_rep <- ifelse(inclusion_table$feature %in% names(M_selected), 1, 0)

                       return(inclusion_rep)
                     }

  end_time <- Sys.time()
  cat("Run time:", end_time - start_time)
  # Stop the parallel backend
  stopCluster(cl)
  
  # Combine the results into the final inclusion table
  for (i in 1:nReps) {
    inclusion_table[, paste0("rep", i)] <- results[, i]
  }
  
  
  ## Calculated the inclusion rate
  column_sums <- ifelse(colSums(inclusion_table[,-1])==0, 1, colSums(inclusion_table[,-1]))
  inclusion_table$inclusion_rate <- rowMeans(sweep(inclusion_table[,-1], 2, column_sums, "/"))
  
  # write.csv(inclusion_table, file = "C:/Users/lclce/Desktop/2nd Sem/Yushu RA/CATSplit/inclusion_table.csv", row.names = FALSE)
  # inclusion_table <- read.csv("C:/Users/lclce/Desktop/2nd Sem/Yushu RA/CATSplit/inclusion_table.csv")
  
  inclusion_rate <- inclusion_table$inclusion_rate
  names(inclusion_rate) <- inclusion_table$feature
  inclusion_rate <- inclusion_rate[order(inclusion_rate, decreasing = FALSE)]
  
  
  ## find the largest l
  cum_sum <- cumsum(inclusion_rate)
  l <- max(which(cum_sum <= 0.05))
  
  ## select features and order them in descending order
  selected_features <- rev(names(inclusion_rate[-(1:l)]))
  
  return(selected_features)
}


CATSplit_noparallel <- function(otu_table, taxonomy_table, meta_table, tree_data, 
                                metric, parallel, nCore, nReps, qval_bound = 0.05){
  inclusion_table <- data.frame(feature = unique(taxonomy_table$taxa_to_genus))
  
  # Non-parallel
  results <- matrix(NA, nrow = nrow(inclusion_table), ncol = nReps)


  for (rep in 1:nReps) {

    ## Data Splitting
    data1_ind <- sample(colnames(otu_table), ncol(otu_table)/2, replace = FALSE)
    data2_ind <- setdiff(colnames(otu_table), data1_ind)

    data1_otu_table <- otu_table[, data1_ind]
    data2_otu_table <- otu_table[, data2_ind]

    data1_meta_table <- meta_table[data1_ind, , drop = FALSE]
    data2_meta_table <- meta_table[data2_ind, , drop = FALSE]

    testList <- unique(taxonomy_table$taxa_to_genus)

    ## Compute R2 values for the first and second halves of the data
    beta1 <- computeR2(testList = testList, otutable = data1_otu_table, taxonomy = taxonomy_table,
                       metaData = data1_meta_table, tree = tree_data, metric = metric,
                       parallel = FALSE, nCore = nCore)[, 3]

    beta2 <- computeR2(testList = testList, otutable = data2_otu_table, taxonomy = taxonomy_table,
                       metaData = data2_meta_table, tree = tree_data, metric = metric,
                       parallel = FALSE, nCore = nCore)[, 3]

    ## Calculate M values and determine tau
    M <- sign(beta1 * beta2) * (abs(beta1) * abs(beta2))
    
    selected_index <- analys(M, abs(M), qval_bound)
    M_selected <- M[selected_index]

    ## Create inclusion indicator vector for this repetition
    inclusion_rep <- ifelse(inclusion_table$feature %in% names(M_selected), 1, 0)

    ## Store the result for this repetition in the results matrix
    results[, rep] <- inclusion_rep
  }
  
  # Combine the results into the final inclusion table
  for (i in 1:nReps) {
    inclusion_table[, paste0("rep", i)] <- results[, i]
  }
  
  
  ## Calculated the inclusion rate
  column_sums <- ifelse(colSums(inclusion_table[,-1])==0, 1, colSums(inclusion_table[,-1]))
  inclusion_table$inclusion_rate <- rowMeans(sweep(inclusion_table[,-1], 2, column_sums, "/"))
  
  # write.csv(inclusion_table, file = "C:/Users/lclce/Desktop/2nd Sem/Yushu RA/CATSplit/inclusion_table.csv", row.names = FALSE)
  # inclusion_table <- read.csv("C:/Users/lclce/Desktop/2nd Sem/Yushu RA/CATSplit/inclusion_table.csv")
  
  inclusion_rate <- inclusion_table$inclusion_rate
  names(inclusion_rate) <- inclusion_table$feature
  inclusion_rate <- inclusion_rate[order(inclusion_rate, decreasing = FALSE)]
  
  
  ## find the largest l
  cum_sum <- cumsum(inclusion_rate)
  l <- max(which(cum_sum <= 0.05))
  
  ## select features and order them in descending order
  selected_features <- rev(names(inclusion_rate[-(1:l)]))
  
  return(selected_features)
}

