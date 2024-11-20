FDR <- NA
Power <- NA
numSelected <- NA
i <- 1

FDR[i] <- mean(sim_res[,1])
Power[i] <- mean(sim_res[,2])
numSelected[i] <- mean(sim_res[,3])
i <- i + 1

df <- data.frame(FDR = FDR, Power = Power, numSelected = numSelected)
rownames(df) <- c("mostAbun_a750_WU", "mostAbun_a100_WU", "leastAbun_a750_WU", "lowCor_a500_WU", "lowCor_a250_WU")

save(df, file = "Perf_df.RData")
