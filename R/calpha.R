calpha_method <- function(casecon, gen) {
  # This function computes C-alpha test  
  nA <- sum(casecon)
  nU <- sum(casecon==0)
  p0 <- nA / (nA + nU)
    
  m <- ncol(gen)
  # copies of the i-th variant type
  n <- apply(gen, 2, function(x) sum(x>0, na.rm=TRUE))
  # copies of the i-th variant type in the cases
  g <- apply(as.matrix(gen[casecon==1,]), 2, function(x) sum(x>0, na.rm=TRUE))
  # Test statistic 
  Talpha <- sum((g - n*p0)^2 - (n * p0 * (1-p0)))
  Talpha
}

# This is C(aplpa) method
calpha <- function(table, permutations=NULL) {  
    casecon <- as.numeric(as.matrix(table[,1]))
    variants <- as.matrix(table[,-1])
    
    ## Checking arguments for correctness:
    if (!is.null(permutations))
    {
      if (mode(permutations) != "numeric" || length(permutations) != 1
          || permutations < 0 || (permutations %% 1) !=0) 
      {
        warning("Argument 'permutations' incorrectly defined. Default value permutations=100 will be used")
        permutations <- 100
      }
    } else {
      permutations <- 0
    } 

    #Total cases
    nA <- length(which(casecon==1))
    #Total controls
    nU <- length(which(casecon==0))
    # % of cases
    p0 <- nA / (nA + nU)
    # Total variants:
    m <- ncol(variants)
    # Copies of the i-th variant type
    n <- apply(as.matrix(variants), 2, function(x) sum(x>0, na.rm=TRUE))
    #print(n)
    # copies of the i-th variant type in the cases
    g <- apply(as.matrix(variants[casecon==1,]), 2, function(x) sum(x>0, na.rm=TRUE))
    
    # Test statistic:
    calpha.stat <- calpha_method(casecon, variants)
    # Variance of Talpha:
    Valpha <- 0
    #for (i in 1:m) {
    #  for (u in 0:n[i]) {
    #    Valpha <- Valpha + (((u - n[i]*p0)^2 - n[i]*p0*(1-p0))^2)*dbinom(u, n[i], p0)
    #  }
    #}
    for (i in 1:m) {
      for (u in 0:n[i]) {
        Valpha <- Valpha + n[i]*(((u - n[i]*p0)^2 - n[i]*p0*(1-p0))^2)*dbinom(u, n[i], p0)
      }
    }
    
    names(Valpha) <- NULL
    # Z-score:
    Zscore <- calpha.stat / sqrt(Valpha)
    
    # p-vaue (asymptotic):
    if (Valpha==0) {
      asym.pval <- 1
    } else {
      asym.pval <- 1 - pchisq(calpha.stat^2 / Valpha, df=1)
    }
    
    # Dealing with permutations:
    p.pval <- NA
    if (permutations != 0)
    {
      x.perm <- rep(0, permutations)
      for (i in 1:permutations)
      {
        perm.sample <- sample(1:length(casecon))
        x.perm[i] <- calpha_method(casecon[perm.sample], variants) 
      }
      # p-value for permutations:
      p.pval <- sum(x.perm^2 > calpha.stat^2) / permutations
    }    
    
    ## Final results
    name <- "C(alpha) Test"
    arg.spec <- c(sum(casecon), 
                  length(casecon)-length(which(casecon==1)), 
                  ncol(variants), 
                  permutations)
    
    names(arg.spec) <- c("cases", "controls", "variants", "n.perms")
    
    res <- list(calpha.stat = calpha.stat, 
               asym.pval = asym.pval, 
               perm.pval = p.pval,
               zscore = Zscore,
               args = arg.spec, 
               name = name)
    
    return(res)
}