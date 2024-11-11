library(glmnet)
library(tidyverse)
library(ape)
library(GUniFrac)
library(MiRKAT)
library(vegan)
library(parallel)
library(doParallel)
source("helper_function.R")
source("preprocess.R")
source("DS_Mann_fun.R")

PoisPara <- 21
sampleSize <- 50
nCore <- 12

meta_table <- data.frame(BinOutcomes = c(rep(1, sampleSize / 2), rep(0, sampleSize / 2)))
rownames(meta_table) <- c(paste("Responder", 1:(sampleSize / 2)), paste("non-Responder", 1:(sampleSize / 2)))

sim_data_file <- paste("input/sim_data_s", sampleSize, "_l", PoisPara, ".RData", sep = "")

load(sim_data_file)

# otu_table <- as.data.frame(sim_data[[1]])
# Mann_WhitU(otu_table, taxonomy_table, meta_table)

spiked_genera <- c("g__Ruminococcus", "g__Veillonella", "g__Lactobacillus",
                   "g__Bacteroides", "g__Prevotella", "g__Holdemania",
                   "g__Escherichia", "g__Olsenella", "g__'Solanum", "g__Collinsella",
                   "g__Bifidobacterium", "g__Corynebacterium", "g__Fusobacterium", "g__Akkermansia", "g__Victivallis", 
                   "g__Methanosphaera", "g__Synergistes", "g__Mucispirillum")
spiked_chain_genera <- unique(taxonomy_table[taxonomy_table$genus %in% spiked_genera,]$taxa_to_genus)

## Mann Whitney
start_time <- Sys.time()
cl <- makeCluster(nCore)
registerDoParallel(cl)
## Parallel Computing
mann_result <- foreach(iter=1:200, .combine = "rbind",
                  .packages=c("ape","GUniFrac","vegan","tidyverse", "MiRKAT")) %dopar% {
                    otu_table <- as.data.frame(sim_data[[iter]])
                    selected_features <- Mann_WhitU(otu_table, taxonomy_table, meta_table, qval_bound = 0.05)
                    
                    fd_pw <- fdp_power(selected_features, spiked_chain_genera[[iter]])
                    num_selected <- length(selected_features)
                    
                    return(c(fd_pw$fdp, fd_pw$power, num_selected))
                  }

stopCluster(cl)
end_time <- Sys.time()
cat("Run time:", end_time - start_time)

colnames(mann_result) <- c("FDR", "Power", "num_selected")


## MDS
start_time <- Sys.time()
cl <- makeCluster(nCore)
registerDoParallel(cl)
## Parallel Computing
DS_result <- foreach(iter=1:200, .combine = "rbind",
                  .packages=c("ape","GUniFrac","vegan","tidyverse", "MiRKAT", "glmnet")) %dopar% {
                    
                    otu_table <- as.data.frame(sim_data[[iter]])
                    selected_features <- DSBin(X = otu_table, y = meta_table$BinOutcomes, taxonomy_table, num_split = 50, qval_bound = 0.05)
                    
                    fd_pw <- fdp_power(selected_features, spiked_chain_genera[[iter]])
                    num_selected <- length(selected_features)
                    
                    return(c(fd_pw$fdp, fd_pw$power, num_selected))
                  }

stopCluster(cl)
end_time <- Sys.time()
cat("Run time:", end_time - start_time)

colnames(DS_result) <- c("FDR", "Power", "num_selected")
