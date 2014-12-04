VT <- function(y, X, maf=0.05, perm=100) {
    ## checking arguments
    Xy_perm = my_check(y, X, perm)
    y = Xy_perm$y
    X = Xy_perm$X
    perm = Xy_perm$perm
    if (mode(maf) != "numeric" || length(maf) != 1
        || maf <= 0  || maf >= 1)
      stop("argument 'maf' must be a value between 0 and 1")
    
    ## running vt method
    mafs = (1 + colSums(X, na.rm=TRUE)) / (2 + 2*nrow(X))
    h.maf = sort(unique(mafs))
    vt.stat = my_vt_method(y, X, mafs, h.maf)
    
    ## permutations
    perm.pval = NA
    if (perm > 0)
    {
      x.perm = rep(0, perm)
      for (i in 1:perm)
      {
        perm.sample = sample(1:length(y))
        x.perm[i] = my_vt_method(y[perm.sample], X, mafs, h.maf)
      }
      ## p-value
      perm.pval = sum(x.perm >= vt.stat) / perm
    }
    
    ## results
    name = "VT: Variable Threshold"
    arg.spec = c(sum(y), length(y)-sum(y), ncol(X), maf, perm)
    names(arg.spec) = c("cases", "controls", "variants", "maf", "n.perms")  
    res = list(vt.stat = vt.stat, 
               perm.pval = perm.pval, 
               args = arg.spec, 
               name = name)
    class(res) = "assoctest"
    return(res)
}

