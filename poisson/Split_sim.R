library(MCMCpack)
library(glmnet)
library(tidyverse)
library(ape)
library(GUniFrac)
library(MiRKAT)
library(vegan)
library(parallel)
library(doParallel)

source("generate_sim.R")
source("Split_fun.R")
source("helper_function.R")
source("preprocess.R")

# define TP genera
spiked_genera <- c("g__Ruminococcus", "g__Veillonella", "g__Lactobacillus",
                   "g__Bacteroides", "g__Prevotella", "g__Holdemania",
                   "g__Escherichia", "g__Olsenella", "g__'Solanum", "g__Collinsella",
                   "g__Bifidobacterium", "g__Corynebacterium", "g__Fusobacterium", "g__Akkermansia", "g__Victivallis",
                   "g__Methanosphaera", "g__Synergistes", "g__Mucispirillum")
spiked_chain_genera <- unique(taxonomy_table[taxonomy_table$genus %in% spiked_genera,]$taxa_to_genus)

# set.seed(123)
# spiked_chain_genera <- list()
# for(i in 1:200){
#   spiked_chain_genera[[i]] <- sample(unique(taxonomy_table$taxa_to_genus), 19, replace = FALSE)
# }
# set.seed(NULL)

## Generate data
marginal<-apply(otutabletoabundance(otu_table),1,mean)
averageCount<-mean(apply(otu_table,2,sum))
testTaxon <- spiked_chain_genera
PoisPara <- 20
sampleSize <- 50
nCore <- 10
metric <- "robust"

meta_table <- data.frame(BinOutcomes = c(rep(1, sampleSize / 2), rep(0, sampleSize / 2)))
rownames(meta_table) <- c(paste("Responder", 1:(sampleSize / 2)), paste("non-Responder", 1:(sampleSize / 2)))

sim_data_file <- paste("input/sim_data_s", sampleSize, "_l", PoisPara, ".RData", sep = "")

cl <- makeCluster(nCore)
registerDoParallel(cl)
if(!file.exists(sim_data_file)){
  sim_data <- foreach(seedNum=1:200,
                      .packages=c("ape","GUniFrac","vegan","dirmult")) %dopar% {
                        generate_sim(seedNum, testTaxon[[seedNum]], PoisPara, taxonomy_table, otu_table, tree_data, averageCount, marginal, sampleSize)
                      }
  save(sim_data, file = sim_data_file)
} else load(sim_data_file)
stopCluster(cl)

## Simulation analysis
start_time <- Sys.time()
cl <- makeCluster(nCore)
registerDoParallel(cl)

sim_res <- foreach(iter=1:50, .combine = "rbind",
                   .packages=c("ape","GUniFrac","vegan","tidyverse", "MiRKAT"),
                   .export = c("compDist", "computeR2", "find_tau")) %dopar% {
                     
                     selected_features <- CATSplit_noparallel(sim_data[[iter]], taxonomy_table, 
                                                    meta_table, tree_data, metric, parallel, nCore, nReps = 50)
                      
                     fd_pw <- fdp_power(selected_features, spiked_chain_genera)
                     num_selected <- length(selected_features)
                     
                     return(c(fd_pw$fdp, fd_pw$power, num_selected))
                   }

stopCluster(cl)
end_time <- Sys.time()
cat("Run time:", end_time - start_time)
save(sim_res, file = "output/sim_res.RData")





##
selected_features <- CATSplit_parallel(sim_data[[1]], taxonomy_table, meta_table, tree_data, metric, parallel, nCore, nReps = 50)
fd_pw <- fdp_power(selected_features, spiked_chain_genera)





# 
# otu_table$taxa_to_genus <- taxonomy_table$taxa_to_genus
# 
# genus_count <- otu_table %>%
#   group_by(taxa_to_genus) %>%
#   summarise(across(where(is.numeric), \(x) sum(x, na.rm = TRUE))) %>%
#   column_to_rownames("taxa_to_genus")
# 
# otu_table <- dplyr::select(otu_table, -taxa_to_genus)
# 
# 
# top20_names <- names(rowMeans(genus_count)[order(rowMeans(genus_count), decreasing = TRUE)])[1:20]
# spiked_table <- genus_count[top20_names,]
# 
# spiked_genera <- sapply(top20_names, function(x) strsplit(x, "\\.")[[1]][6])
# names(spiked_genera) <- NULL
# 
# # plot selected genera distance
# selected_genera <- spiked_genera
# spiked_genera %in% taxonomy_table$genus
# 
# taxonomy <- taxonomy_table
# 
# 
# selected_OTUs <- rownames(taxonomy[taxonomy$genus %in% spiked_genera,])
# 
# # 2. Subset the Phylogenetic Tree for the Selected OTUs
# selected_tree <- drop.tip(tree_data, setdiff(tree_data$tip.label, selected_OTUs))
# 
# # 3. Calculate Pairwise Phylogenetic Distances
# dist_matrix <- cophenetic(selected_tree)
# 
# # 4. Map Each OTU to Its Genus
# taxonomy$OTU_ID <- rownames(taxonomy)  # Add OTU IDs if not present
# selected_taxonomy <- taxonomy[taxonomy$OTU_ID %in% selected_tree$tip.label,]
# genus_map <- selected_taxonomy[, c("OTU_ID", "genus")]
# 
# # 5. Aggregate Distances at the Genus Level
# genus_distance <- matrix(NA, nrow = length(selected_genera), ncol = length(selected_genera))
# rownames(genus_distance) <- colnames(genus_distance) <- selected_genera
# 
# for (genus1 in selected_genera) {
#   for (genus2 in selected_genera) {
#     otus1 <- genus_map[genus_map$genus == genus1, "OTU_ID"]
#     otus2 <- genus_map[genus_map$genus == genus2, "OTU_ID"]
#     if (length(otus1) > 0 & length(otus2) > 0) {
#       # Compute the mean pairwise distance between OTUs of the two genera
#       genus_distance[genus1, genus2] <- mean(dist_matrix[otus1, otus2], na.rm = TRUE)
#     }
#   }
# }
# 
# # 6. Visualize the Genus Distance Matrix as a Heatmap
# heatmap(as.matrix(genus_distance), symm = TRUE, col = heat.colors(256), main = "Genus-Level Phylogenetic Distances")
# 
# # 7. Optional: Visualize with a Dendrogram for Clustering of Genera
# dist_genus <- as.dist(genus_distance)
# hclust_genus <- hclust(dist_genus, method = "average")
# plot(hclust_genus, main = "Dendrogram of Selected Genera", xlab = "Genera", ylab = "Phylogenetic Distance")

