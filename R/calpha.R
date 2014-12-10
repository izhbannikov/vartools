calpha_method <- function(casecon, gen) {
    # Internal function for BST method
    nA = sum(casecon)
    nU = sum(casecon==0)
    p0 = nA / (nA + nU)
    
    m = ncol(gen)
    # copies of the i-th variant type
    n = apply(gen, 2, function(x) sum(x>0, na.rm=TRUE))
    # copies of the i-th variant type in the cases
    g = apply(gen[casecon==1,], 2, function(x) sum(x>0, na.rm=TRUE))
    # Test statistic 
    Talpha = sum((g - n*p0)^2 - (n * p0 * (1-p0)))
    Talpha
}

# This is C(aplpa) method
calpha <- function(table, perm=NULL) {  
    y <- as.numeric(as.matrix(table[,1]))
    X <- as.matrix(table[,-1])
    ## checking arguments
    if (!is.null(perm))
    {
      if (mode(perm) != "numeric" || length(perm) != 1
          || perm < 0 || (perm %% 1) !=0) 
      {
        warning("Argument 'perm' incorrectly defined. Value perm=100 is used")
        perm = 100
      }
    } else perm=0
    Xy_perm = check_args(y, X, perm)
    y = Xy_perm$y
    X = Xy_perm$X
    perm = Xy_perm$perm
    
    # how many cases
    nA = sum(y)
    # how many controls
    nU = sum(y == 0)
    # proportion of cases
    p0 = nA / (nA + nU)
    # how many variants
    m = ncol(X)
    # copies of the i-th variant type
    n = apply(X, 2, function(x) sum(x>0, na.rm=TRUE))
    # copies of the i-th variant type in the cases
    g = apply(X[y==1,], 2, function(x) sum(x>0, na.rm=TRUE))
    
    # Test statistic 
    calpha.stat = calpha_method(y, X)
    # Variance of Talpha
    Valpha = 0
    for (i in 1:m) {
      for (u in 0:n[i]) {
        Valpha = Valpha + (((u - n[i]*p0)^2 - n[i]*p0*(1-p0))^2)*dbinom(u, n[i], p0)
      }
    }
    names(Valpha) = NULL
    # Z score
    Zscore = calpha.stat / sqrt(Valpha)
    
    # asymptotic p-vaue
    if (Valpha==0) asym.pval=1 else
      asym.pval = 1 - pchisq(calpha.stat^2 / Valpha, df=1)
    
    # permutations
    perm.pval = NA
    if (perm != 0)
    {
      x.perm = rep(0, perm)
      for (i in 1:perm)
      {
        perm.sample = sample(1:length(y))
        x.perm[i] = calpha_method(y[perm.sample], X) 
      }
      # p-value 
      perm.pval = sum(x.perm^2 > calpha.stat^2) / perm
    }    
    
    ## results
    name = "C(alpha) Test"
    arg.spec = c(sum(y), length(y)-sum(y), ncol(X), perm)
    names(arg.spec) = c("cases", "controls", "variants", "n.perms")
    res = list(calpha.stat = calpha.stat, 
               asym.pval = asym.pval, 
               perm.pval = perm.pval,
               zscore = Zscore,
               args = arg.spec, 
               name = name)
    return(res)
}

