rarecov_method <- function(casecon, gen, dif) {
  #casecon <- y
  #gen <- X.new
  #dif <- 0.5
  # Internal function for RARECOV method
  # rare cover algorithm
  M <- ncol(gen)
  selected <- NULL
  temp.stats <- rep(0, M)
  temp.pvals <- rep(0, M)
    
  ## get the first variant
  for (j in 1:M) {
    temp <- chisq.test(table(gen[,j], casecon),simulate.p.value = T, B=2)
    temp.stats[j] <- temp$statistic
    temp.pvals[j] <- temp$p.value
  }
  win.j <- which(temp.stats == max(temp.stats))
  selected <- c(selected, win.j[1])
  Xcorr <- temp.stats[win.j[1]]
  Xpval <- temp.pvals[win.j[1]]
  rest <- setdiff(1:M, selected)
    
  ## get rest of variants
  Xcorr.dif <- 1
  while (Xcorr.dif > dif){   
    temp.stats <- rep(0, M) 
    temp.pvals <- rep(0, M)
    for (j in 1:length(rest)){
        gen.new <- rowSums(as.matrix(gen[,c(selected,rest[j])]), na.rm=TRUE)
        gen.new[gen.new != 0] <- 1
        temp <- chisq.test(table(gen.new, casecon),simulate.p.value=T)
        temp.stats[rest[j]] <- temp$statistic
        temp.pvals[rest[j]] <- temp$p.value
    }
    win.j <- which(temp.stats == max(temp.stats))
    Xcorr.new <- temp.stats[win.j[1]]
    Xpval.new <- temp.pvals[win.j[1]]
    Xcorr.dif <- Xcorr.new - Xcorr
    if (Xcorr.dif > dif){
        selected <- c(selected, win.j)
        Xcorr <- Xcorr.new
        Xpval <- Xpval.new
        rest <- setdiff(1:M, selected)
    }
  }
  list(stat=Xcorr, pval=Xpval, sel=selected)
}

rarecover <- function(table, maf=0.05, dif=0.5, perm=100) {
  y <- as.numeric(as.matrix(table[,1]))
  X <- as.matrix(table[,-1])
  ## checking arguments
  Xy_perm <- check_args(y, X, perm)
  y <- Xy_perm$y
  X <- Xy_perm$X
  perm <- Xy_perm$perm
  
  if (mode(maf)!= "numeric" || length(maf) != 1 || maf<=0 || maf>1)
      stop("argument 'maf' incorreclty defined; must be a value between 0 and 1")
  if (mode(dif) != "numeric" || length(dif) != 1
        || dif <= 0  || dif >= 1)
      stop("argument 'dif' must be a value between 0 and 1")
    
  ## get minor allele frequencies
  MAF <- colMeans(X, na.rm=TRUE) / 2   
  ## how many variants < maf
  rare <- sum(MAF < maf)
  print(MAF)
  print(rare)
  if (rare == 0)
      stop(paste("Error:  no rare variants detected below maf=", maf, sep=""))
  ## collapsing
  if (rare == 1){   
    # if rare variants <= 1, then NO collapse is needed
    X.new <- X
  } else {
    X.new <- X[,MAF < maf]     
  }
  ## change genotype 2 into 1
  X.new[X.new == 2] <- 1
  rarecov <- rarecov_method(y, X.new, dif)
  rc.stat <- rarecov$stat
  chisq.pval <- rarecov$pval
  rc.sel <- rarecov$sel
  names(rc.sel) <- NULL
    
  ## permutations
  perm.pval <- NA
  if (perm > 0){
    x.perm <- rep(0, perm)
    for (i in 1:perm){
      perm.sample <- sample(1:length(y))
      rarecov.perm <- rarecov_method(y[perm.sample], X.new, dif)
      x.perm[i] <- rarecov.perm$stat
    }
    ## p-value
    perm.pval <- sum(x.perm > rc.stat) / perm
  }
    
  ## results
  name <- "RareCover Test"
  arg.spec <- c(sum(y), length(y)-sum(y), ncol(X), rare, maf, length(rc.sel), perm)
  arg.spec <- as.character(arg.spec)
  names(arg.spec) <- c("cases", "controls", "variants", "rarevar", "maf", "varsel", "n.perms")
  sel.names <- names(X)[rc.sel] 
  if (is.null(sel.names)) 
      sel.names <- paste("var", rarecov$sel, sep="")
    names(rc.sel) <- sel.names
    res <- list(rc.stat = rc.stat, 
               perm.pval = perm.pval, 
               set = rc.sel, 
               args = arg.spec, 
               name = name)
    return(res)
}