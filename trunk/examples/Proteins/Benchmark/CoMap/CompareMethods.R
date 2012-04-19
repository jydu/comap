# Created on 17/04/2012 by jdutheil
# Compare mapping methods


compare.mapping<-function(map1, map2) {
  layout(matrix(1:2, nrow=1))
  x<-unlist(map1[, 3:ncol(map1)])
  y<-unlist(map2[, 3:ncol(map2)])
  plot(y~x, pch=19)
  abline(0,1, col="red")
  hist(y-x)
  abline(v=0, col="red")
}

noweight.naive<-read.table("Myo_naive.vec", header=TRUE, sep="\t")
noweight.lapla<-read.table("Myo_laplace.vec", header=TRUE, sep="\t")
noweight.unifo<-read.table("Myo_unif.vec", header=TRUE, sep="\t")
noweight.decom<-read.table("Myo_decomp.vec", header=TRUE, sep="\t")

compare.mapping(noweight.decom, noweight.lapla)
compare.mapping(noweight.decom, noweight.unifo)
compare.mapping(noweight.decom, noweight.naive)

weight.naive<-read.table("Myo_naive_grantham.vec", header=TRUE, sep="\t")
weight.unifo<-read.table("Myo_unif_grantham.vec", header=TRUE, sep="\t")
weight.decom<-read.table("Myo_decomp_grantham.vec", header=TRUE, sep="\t")

compare.mapping(weight.decom, weight.unifo)
compare.mapping(weight.decom, weight.naive)

compare.mapping(weight.decom, noweight.decom)
compare.mapping(weight.unifo, noweight.unifo)
compare.mapping(weight.naive, noweight.naive)
