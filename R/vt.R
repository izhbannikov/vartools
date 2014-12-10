vt_method <- function(casecon, gen, mafs, h.maf) {
    # mafs: minor allele frequencies
    # h.maf: unique mafs
    z.scores = rep(0, length(h.maf)-1)
    y.new = casecon - mean(casecon)
    for (i in 1:(length(h.maf)-1))
    {
      z.num = sum(gen[,mafs<h.maf[i+1]] * y.new, na.rm=TRUE)
      z.denom = sqrt(sum((gen[,mafs<h.maf[i+1]])^2, na.rm=TRUE))
      z.scores[i] = z.num / z.denom
    }
    stat = max(z.scores)
    stat
}

vt <- function(table, maf=0.05, perm=50) {
    y <- as.numeric(as.matrix(table[,1]))
    X <- as.matrix(table[,-1])
    ## checking arguments
    Xy_perm = check_args(y, X, perm)
    y = Xy_perm$y
    X = Xy_perm$X
    perm = Xy_perm$perm
    if (mode(maf) != "numeric" || length(maf) != 1
        || maf <= 0  || maf >= 1)
      stop("argument 'maf' must be a value between 0 and 1")
    
    ## running vt method
    mafs = (1 + colSums(X, na.rm=TRUE)) / (2 + 2*nrow(X))
    h.maf = sort(unique(mafs))
    vt.stat = vt_method(y, X, mafs, h.maf)
    
    ## permutations
    perm.pval = NA
    if (perm > 0)
    {
      x.perm = rep(0, perm)
      for (i in 1:perm)
      {
        perm.sample = sample(1:length(y))
        x.perm[i] = vt_method(y[perm.sample], X, mafs, h.maf)
      }
      ## p-value
      perm.pval = sum(x.perm >= vt.stat) / perm
    }
    
    ## results
    name = "Variable Threshold (VT) method"
    arg.spec = c(sum(y), length(y)-sum(y), ncol(X), maf, perm)
    names(arg.spec) = c("cases", "controls", "variants", "maf", "n.perms")  
    res = list(vt.stat = vt.stat, 
               perm.pval = perm.pval, 
               args = arg.spec, 
               name = name)
    return(res)
}

