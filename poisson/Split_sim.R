library(MCMCpack)
library(glmnet)
library(tidyverse)
library(ape)
library(GUniFrac)
library(MiRKAT)
library(vegan)
library(parallel)
library(doParallel)
library(stringr)
library(dirmult)

source("generateSim1.R")
source("Split_fun.R")
source("helper_function.R")
source("preprocess.R")
load("../mostLeastAbun.RData")
#load("mostVarLeastVar.RData")
#load("highLowCor.RData")
# for mostVar, find the value after the last .
#spiked_chain_genera <- leastVar
#spiked_chain_genera<-mostVar
spiked_chain_genera <- mostAbun # 9.1% FDR 100% power with euclidean
#spiked_chain_genera <- leastAbun 100 percent 
#spiked_chain_genera <- highCor
#spiked_chain_genera <- lowCor
str(otu_table)
str(taxonomy_table)

identical(rownames(otu_table), rownames(taxonomy_table))
## Generate data
marginal<-apply(otutabletoabundance(otu_table),1,mean)
averageCount<-mean(apply(otu_table,2,sum))
testTaxon <- spiked_chain_genera
addedCount <- 100
sampleSize <- 44
#metric <- "Weighted UniFrac"
metric <- "euclidean"
DiriSum<-62

# count how many OTU under the leastVar taxa
for(i in 1:length(spiked_chain_genera)){
  print(sum(taxonomy_table$taxa_to_genus == spiked_chain_genera[i]))
}

meta_table <- data.frame(BinOutcomes = c(rep(1000, sampleSize / 2), rep(0, sampleSize / 2)))
rownames(meta_table) <- c(paste("Responder", 1:(sampleSize / 2)), paste("non-Responder", 1:(sampleSize / 2)))

OTUsim<-generate_sim(1, testTaxon, addedCount, taxonomy_table, otu_table, tree_data, averageCount, marginal, sampleSize,DiriSum)
otu_table<-OTUsim
inclusion_table <- data.frame(feature = unique(taxonomy_table$taxa_to_genus))

#selected_features <- singleSplit(inclusion_table, OTUsim, taxonomy_table, 
#                                       meta_table, tree_data, metric)
#selected_features<-inclusion_table$feature[which(selected_features>0)]

selected_features <- CATSplit_parallel(OTUsim, taxonomy_table, 
                                        meta_table, tree_data, metric, nCore=24, nReps = 50)
length(selected_features)
1-sum(selected_features%in% spiked_chain_genera)/length(selected_features) # FDP
sum(spiked_chain_genera%in% selected_features)/length(spiked_chain_genera) # Power)                      

# fd_pw <- fdp_power(selected_features, spiked_chain_genera)
# num_selected <- length(selected_features)
                     
