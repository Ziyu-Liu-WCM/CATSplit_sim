#otutable to otu abundance

otutabletoabundance<-function(otutable){
colsum<-apply(otutable,2,sum)
colmax<-matrix(rep(colsum,each=nrow(otutable)),nrow=nrow(otutable))
otu.abundance.table<-otutable/colmax
}
