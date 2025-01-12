library(foreach)
library(doParallel)
setwd("~/Desktop/temp")
permR2<-function(X,j,y){
  Xperm<-X
  Xperm[,j]<-sample(X[,j],replace=FALSE)
  data<-data.frame(y=y,Xperm)
  fit_perm<-lm(y~.,data=data)
  return(summary(fit_perm)$r.squared)
}

testFun<-function(i){
  set.seed(i)
  n <-800
  p <- 250
  p0 <- 20
  X <- mvrnorm(n, mu = rep(0, p), Sigma = diag(p))
  y<-rnorm(n)
  X=as.data.frame(X)
  names(X)=paste0('X',1:p)
  
  sample_index1 <- sample(x = c(1:n), size = 0.5 * n, replace = F)
  sample_index2 <- setdiff(c(1:n), sample_index1)
  lm1<-lm(y[sample_index1]~as.matrix(X[sample_index1,]))
  R2orig1<-summary(lm1)$r.squared
  lm2<-lm(y[sample_index2]~as.matrix(X[sample_index2,]))
  R2orig2<-summary(lm2)$r.squared
  beta1<-sapply(1:ncol(X),function(j) permR2(X[sample_index1,],j,y[sample_index1]))-R2orig1
  beta2<-sapply(1:ncol(X),function(j) permR2(X[sample_index2,],j,y[sample_index2]))-R2orig2
  
  mirror<-sign(beta1*beta2)*(abs(beta1)+abs(beta2))
  return(mirror)
}

  cl <- makeCluster(detectCores(), outfile='LOG.TXT')
  registerDoParallel(cl)
  MirrorStat<-foreach(seedNum=201:400,
          .packages=c("MASS")) %dopar% {
            testFun(seedNum)
          }

  stopImplicitCluster()
  stopCluster(cl)
  # Deregister the parallel backend
  registerDoSEQ()

  # Remove objects
  rm(cl)
  #Call garbage collection
  gc()

# CONVERT A LIST TO A MATRIX FOR MIRRORSTAT
MirrorStat<-do.call(rbind,MirrorStat)
str(MirrorStat)
sum(apply(MirrorStat,1,median)>0)
