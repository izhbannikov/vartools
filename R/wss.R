wss_method <- function(casecon, gen) {
  controls <- !casecon
  n.loc <- ncol(gen)
  
  ## calculate weights
  w <- rep(0, n.loc)
  for (j in 1:n.loc) {
    nNA <- is.na(gen[,j])
    mU <- sum(gen[controls,j], na.rm=TRUE)
    nU <- sum(controls[!nNA])  	
    q <- (mU+1) / (2*nU+2)
    n <- sum(!nNA)
    w[j] <- sqrt(n * q * (1-q))
  }
  score <- NULL
  ## calculate genetic score
  if (length(w) > 1) {
    score <- rowSums(gen %*% diag(1/w), na.rm=TRUE) 
  } else {
    score <- rowSums(gen %*% w, na.rm=TRUE) 
  }
  # rank.score = order(score)
  rank.score <- rank(score)
  
  ## sum of ranks of cases
  x <- sum(rank.score[which(casecon==1)])
  
  return(x)
}

wss <- function(table, perm=100) {
  
  y <- as.numeric(as.matrix(table[,1]))
  X <- as.matrix(table[,-1])
  
  
  ## checking arguments
  Xy_perm <- check_args(y, X, perm)
  y <- Xy_perm$y
  X <- Xy_perm$X
  perm <- Xy_perm$perm
  
  ## number of individuals N
  N = nrow(X)
  ## number of variants
  M = ncol(X)
    
  ## running wss method
  wss.stat <- wss_method(y, X)  
  
    
  ## permutations
  perm.pval <- NA
  asym.pval <- NA
  if (perm > 0) {
    x.perm <- rep(0, perm)
    for (i in 1:perm) {
      perm.sample <- sample(1:length(y))
      x.perm[i] <- wss_method(y[perm.sample], X) 
    }
    # p-value 
    perm.pval <- sum(x.perm > wss.stat) / perm
  } else {
    ## Asymptotic p-values
    # under the null hypothesis test statistics follows a standard normal distribution 
    x.perm <- rep(0, N)
    for (i in 1:N) {
      perm.sample <- sample(1:length(y))
      x.perm[i] <- wss_method(y[perm.sample], X) 
    }
    # p-value 
    z.score <- (wss.stat - mean(x.perm))/sd(x.perm)
    asym.pval = 1 - pnorm(z.score)
  }
    
  ## results
  name <- "WSS: Weighted Sum Statistics"
  arg.spec <- c(sum(y), length(y)-sum(y), ncol(X), perm)
  names(arg.spec) <- c("cases", "controls", "variants", "n.perms")  
  res <- list(wss.stat = wss.stat, 
              asym.pval = asym.pval,
              perm.pval = perm.pval, 
              args = arg.spec, 
              name = name)
  return(res)
}