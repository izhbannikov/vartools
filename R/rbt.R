rbt_method <- function(casecon, gen)
  {
    # Internal function for RBT method
    ## num of counts of rare allele in cases A and controls U
    fA = colSums(gen[casecon==1,], na.rm=TRUE)
    fU = colSums(gen[casecon==0,], na.rm=TRUE)
    
    ## RBT S+ (enrichment of mutations in cases)
    ## get those A larger than U, and sort them
    AltU = fA > fU  # get those A larger than U
    if (sum(AltU) > 0)
    {
      kA.plus = sort(fA[AltU])
      kU.plus = fU[AltU][order(fA[AltU])]
      ## how many unique A > U
      uniqA.plus = unique(sort(fA[AltU]))
      uniqU.plus = unique(sort(fU[AltU]))
      if (length(uniqA.plus) == 1)
      {
        fua = (uniqU.plus + uniqA.plus) / 2
        pU = ppois(uniqU.plus, fua)
        pA = 1 - ppois(uniqA.plus - 1, fua) 
        PUA.plus = -log(pU * pA)
        NUA.plus = 1
      } else {
        ## matrix of observed counts
        ncolsA.plus = length(uniqA.plus)
        nrowsU.plus = length(uniqU.plus)  	
        NUA.plus = matrix(0, nrowsU.plus, ncolsA.plus)
        
        for (j in 1:ncolsA.plus) {
          for (i in kU.plus[kA.plus==uniqA.plus[j]]) {
            NUA.plus[which(uniqU.plus==i),j] = sum(kU.plus[kA.plus==uniqA.plus[j]] == i)
          }
        }
        # dimnames(NUA.plus) = list(uniqU.plus, uniqA.plus)
        
        # matrix of poisson probabilities
        PUA.plus = matrix(0, nrowsU.plus, ncolsA.plus)
        for (i in 1:nrowsU.plus) {
          for (j in 1:ncolsA.plus) {
            fua = (uniqU.plus[i] + uniqA.plus[j]) / 2
            pU = ppois(uniqU.plus[i], fua)
            pA = 1 - ppois(uniqA.plus[j]-1, fua) 
            PUA.plus[i,j] = -log(pU * pA)
          }
        }
        # dimnames(PUA.plus) = list(uniqU.plus, uniqA.plus)
      }
    } else { # sum(AltU) <= 0
      NUA.plus = 0
      PUA.plus = 0
    } 
    
    ## RBT S- (enrichment of mutations in controls)
    ## get those A smaller than U, and sort them
    AstU = fA < fU  # get those A smaller than U
    if (sum(AstU) > 0)
    {
      kA.minus = sort(fA[AstU])
      kU.minus = fU[AstU][order(fA[AstU])]
      ## how many unique A > U
      uniqA.minus = unique(sort(fA[AstU]))
      uniqU.minus = unique(sort(fU[AstU]))
      if (length(uniqA.minus) == 1)
      {
        fua = (uniqU.minus + uniqA.minus) / 2
        pU = ppois(uniqU.minus, fua)
        pA = 1 - ppois(uniqA.minus - 1, fua) 
        PUA.minus = -log(pU * pA)
        NUA.minus = 1
      } else {				
        ## matrix of observed counts
        ncolsA.minus = length(uniqA.minus)
        nrowsU.minus = length(uniqU.minus)
        NUA.minus = matrix(0, nrowsU.minus, ncolsA.minus)
        
        for (j in 1:ncolsA.minus) {
          for (i in kU.minus[kA.minus==uniqA.minus[j]]) {
            NUA.minus[which(uniqU.minus==i),j] = sum(kU.minus[kA.minus==uniqA.minus[j]] == i)
          }
        }
        # dimnames(NUA.minus) = list(uniqU.minus, uniqA.minus)
        
        # matrix of poisson probabilities
        PUA.minus = matrix(0, nrowsU.minus, ncolsA.minus)
        for (i in 1:nrowsU.minus) {
          for (j in 1:ncolsA.minus) {
            fua = (uniqU.minus[i] + uniqA.minus[j]) / 2
            pU = ppois(uniqU.minus[i], fua)
            pA = 1 - ppois(uniqA.minus[j]-1, fua) 
            PUA.minus[i,j] = -log(pU * pA)
          }
        }
        # dimnames(PUA.minus) = list(uniqU.minus, uniqA.minus)
      }
    } else {
      NUA.minus = 0
      PUA.minus = 0
    }
    
    ## RBT S+ (enrichment of mutations in cases)
    rbt.plus = sum(NUA.plus * PUA.plus)
    ## RBT S- (enrichment of mutations in controls)
    rbt.minus = sum(NUA.minus * PUA.minus)
    ## RBT statistic	
    stat = max(rbt.plus, rbt.minus)
    stat
}

rbt <- function(table, perm=150) {
    y <- as.numeric(as.matrix(table[,1]))
    X <- as.matrix(table[,-1])
    ## checking arguments
    Xy_perm = check_args(y, X, perm)
    y = Xy_perm$y
    X = Xy_perm$X
    perm = Xy_perm$perm
    
    ## apply RBT method
    rbt.stat = rbt_method(y, X)
    
    ## permutations
    perm.pval = NA
    if (perm > 0)
    {
      x.perm = rep(0, perm)
      for (i in 1:perm)
      {
        perm.sample = sample(1:length(y))
        x.perm[i] = rbt_method(y[perm.sample], X) 
      }
      # p-value 
      perm.pval = sum(x.perm > rbt.stat) / perm
    }
    
    ## results
    name = "Replication Based Test (RBT) method"
    arg.spec = c(sum(y), length(y)-sum(y), ncol(X), perm)
    names(arg.spec) = c("cases", "controls", "variants", "n.perms")	
    res = list(rbt.stat = rbt.stat, 
               perm.pval = perm.pval, 
               args = arg.spec, 
               name = name)
    return(res)
}

