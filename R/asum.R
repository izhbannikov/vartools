asum_method <- function(U, V) {
  # Internal function for ASUM method
  # decreasing order of U
  Ords <- order(U/sqrt(diag(V)), decreasing=TRUE)
  # Ustd in decreasing order
  Uords <- U[Ords]
  # ord cov matrix
  V.ords <- as.matrix(V[Ords, Ords])
  
  ## get scores
  k <- length(U) 
  scores <- scores.ord <- rep(0, k)
  pvals <- pvals.ord <- rep(0, k)
  for (j in 1:k){
      aux <- sum_method(U[1:j], V[1:j,1:j])
      scores[j] <- aux[1]
      pvals[j] <- aux[2]
      aux.ord <- sum_method(Uords[1:j], V.ords[1:j,1:j])
      scores.ord[j] <- aux.ord[1]
      pvals.ord[j] <- aux.ord[2]
  }
  p1 <- min(pvals)
  p2 <- min(pvals.ord)
  S1 <- scores[which(pvals==p1)]
  S2 <- scores.ord[which(pvals.ord==p2)] 
  
  # Return value
  c(S1, p1, S2, p2)
}

asum <- function(table, perm=100) {
  y <- as.numeric(as.matrix(table[,1]))
  X <- as.matrix(table[,-1])  
  ## checking arguments
  Xy_perm <- check_args(y, X, perm)
  y <- Xy_perm$y
  X <- Xy_perm$X
  perm <- Xy_perm$perm
    
  ## get U and V
  getuv <- getUV(y, as.matrix(X))
  U <- getuv$U
  V <- getuv$V
  ## run score method
  stat.asum <- asum_method(U, V)
  asum.stat1 <- stat.asum[1]  # stat normal
  p1.asum <- stat.asum[2]     # pval normal
    
  ## permutations
  perm.pval <- NA
  if (perm > 0){
      p1.perm <- rep(0, perm)
      ymean <- mean(y)
      for (i in 1:perm) {
        perm.sample <- sample(1:length(y))
        # center phenotype y
        y.perm <- y[perm.sample] - ymean
        # get score vector
        U.perm <- colSums(y.perm * X, na.rm=TRUE)
        perm.asum <- asum_method(U.perm, V)
        p1.perm[i] <- perm.asum[2]
      }
      # p-value 
      perm.pval <- sum(p1.perm > p1.asum) / perm   # normal
  }
    
   ## results
    name <- "aSum: Adaptive Sum Test"
    arg.spec <- c(sum(y), length(y)-sum(y), ncol(X), perm)
    names(arg.spec) <- c("cases", "controls", "variants", "n.perms")
    res <- list(asum.stat = asum.stat1, 
               perm.pval = perm.pval, 
               args = arg.spec, 
               name = name)
    return(res)
}