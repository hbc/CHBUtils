
qval_asymm <- function(ts, ps){
  require(qvalue)
  m <- length(ps)  ##total number of genes
  m1 <- length(ps[ts <= 0]) ##number of genes with negative test statistic
  m2 <- length(ps[ts > 0])  ##number of genes with positive test statistic
  p1 <- ps[ts <= 0]  ##p-values for genes with negative test statistics
  p2 <- ps[ts > 0]  ##p-values for genes with positive test statistics

  p0 <- qvalue(ps)$pi0  ##estimate of pi0, proportion of genes that are EE
  m0 <- p0*m  ##estimate of number of EE genes
  m0half <- m0/2 ##estimate of number of EE with negative (positive) signs
   
  rkp1 <- rank(p1)
  rkp2 <- rank(p2)

  fdrh1 <- p1*m0half/rkp1
  fdrh2 <- p2*m0half/rkp2

  qv1 <- sapply(rkp1, function(x) min(fdrh1[rkp1>=x]))
  qv2 <- sapply(rkp2, function(x) min(fdrh2[rkp2>=x]))

  qvasymm <- rep(NA, m)
  qvasymm[ts <= 0] <- qv1
  qvasymm[ts > 0] <- qv2

  qvasymm
}

