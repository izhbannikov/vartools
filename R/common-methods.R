sum_method <- function(U, V) {
    # score statistic and p-value
    if (abs(sum(U)) < 1e-20) {
      score = 0
      pval = 1
    } else {
      ones = rep(1, length(U))
      score = sum(U) / sqrt(sum(V))
      pval = 1 - pchisq(score^2, 1)
    }
    res = c(score, pval)
    res
}

getUV <- function(y, X) {
    # Internal function for getUV
    # center phenotype y
    y.new = y - mean(y)
    # get score vector
    U = colSums(y.new * X, na.rm=TRUE)
    
    # get covariance matrix
    X.new = scale(X, scale=FALSE)
    if (sum(complete.cases(X)) != length(y))  { # missing data
      tX.new = t(X.new)
      Sx = matrix(0, ncol(X), ncol(X))
      for (i in 1:nrow(tX.new))
      {
        for (j in i:ncol(X.new))
        {
          Sx[i,j] = sum(tX.new[i,] * X.new[,j], na.rm=TRUE)
        }
      }
      sx.diag = diag(Sx)
      Sx = Sx + t(Sx)
      diag(Sx) = sx.diag
      V = mean(y) * (1 - mean(y)) * Sx
    } else {  # no missing data
      V = mean(y) * (1 - mean(y)) * (t(X.new) %*% X.new)
    }
    # results
    res.uv = list(U=U, V=V)
    return(res.uv)
}


check_args <- function(y, X, perm) {
    ## Internal function to check argument:s y, X, perm
    
    # y as numeric vector
    if (!is.vector(y) || mode(y) != "numeric")
      stop("Sorry, argument 'y' must be a numeric vector")
    # no misssing values in y
    if (any(is.na(y))) 
      stop("No missing data allowed in argument 'y' ")  
    # binary values (0, 1) in y
    if (!all(y %in% c(0, 1)))
      stop("Argument 'y' must contain only 0 and 1")
    # X as matrix or data frame
    if(!is.matrix(X) & !is.data.frame(X))
      stop("Argument 'X' must be a matrix or data.frame")    
    # compatibility between X and y
    if (nrow(X) != length(y)) 
      stop("'X' and 'y' have different lengths")
    # force X as matrix
    if (!is.matrix(X)) X = as.matrix(X)
    # permutations
    if (mode(perm) != "numeric" || length(perm) != 1
        || perm < 0 || (perm %% 1) !=0) 
    {
      warning("argument 'perm' incorrectly defined. Value perm=100 is used")
      perm = 100
    }
    # results
    list(y=y, X=X, perm=perm)
  }