generate_sim <- function(iter, testTaxon, FoldChange, taxonomy_table, otu_table, tree, averageCount, marginal, sampleSize, scalefac){
  
  set.seed(iter)
  
  OTUsim <- matrix(0, ncol = sampleSize, nrow = nrow(otu_table))
  rownames(OTUsim) <- rownames(otu_table)
  
  i <- 1
  for (i in 1:ncol(OTUsim)) {
    OTUsim[, i] <- rmultinom(1, averageCount, rdirichlet(1, marginal * scalefac))
    if (FoldChange > 1) {
      if (i <= (sampleSize / 2)) {
        OTUsim[taxonomy_table$taxa_to_genus %in% testTaxon, i] <- OTUsim[taxonomy_table$taxa_to_genus %in% testTaxon, i] * FoldChange
      }
    }
  }
  
  
  colnames(OTUsim) <- c(paste("Responder", 1:(sampleSize / 2)), paste("non-Responder", 1:(sampleSize / 2)))
  
  return(OTUsim)
}


