kbac <- function(table, alpha=NULL, num.permutation=100, quiet = T, maf.upper.bound = 1.0, alternative = 1) {
    ydatIn <- as.matrix(table[,1])
    xmat <- as.matrix(table[,-1])
    mafIn <- apply(xmat, 2, function(x) sum(x[which(x > 0)])) / (length(ydatIn) * 2)
    xmat[which(xmat != 1 & xmat != 2 & xmat != 0)] <- 0
    xdatIn <- matrix(t(xmat), nrow = 1)
    xcol <- ncol(xmat)
    ylen <- nrow(xmat)
    nn <- num.permutation
    aa <- alpha
    qq <- quiet
    mafUpper <- maf.upper.bound
    test_results <- KbacGetP(nn, qq, aa, mafUpper, xdatIn, ydatIn, mafIn, xcol, ylen, alternative)
    
    ## results
    name <- "KBAC: kernel-based adaptive cluster"
    arg.spec <- c(sum(table[,1]), length(table[,1])-sum(table[,1]), xcol, num.permutation)
    names(arg.spec) <- c("cases", "controls", "variants", "n.perms")  
    res <- list(kbac.stat = test_results$kbac.stat, 
                asym.pval = ifelse(num.permutation==0, test_results$pvalue, NA),
                perm.pval = ifelse(num.permutation==0, NA, test_results$pvalue), 
                args = arg.spec, 
                name = name)
    
    return(res)
}
