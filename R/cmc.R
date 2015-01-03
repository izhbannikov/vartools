cmc_method <- function(casecon, gen, method) {
  markers <- casecon
  status <- ((rowSums(as.matrix(gen)))>0)
  
  # Contingency table:
  ct <- matrix(nrow=2,ncol=2)
  ct[1,1] <- sum(markers>=1 & status==1)
  ct[1,2] <- sum(markers==0 & status==1)
  ct[2,1] <- sum(markers>=1 & status==0)
  ct[2,2] <- sum(markers==0 & status==0)
  
  if (method == "fisher") {
    stat <- fisher.test(ct) #$p.value # Fisher's exact test
  }
  if (method == "chisq") {
    stat <- chisq.test(ct) #$p.value # Should this test work?
  }
  
  stat
  
}

cmc <- function(table, method="fisher", perm=100, print=F) {
  #cmc <- function(table, maf=0.05, method="fisher", perm=100, print=F) {
  # Input table consists of columns and rows. Rows are individuals and columns are:
  # First column represents a status (disease / no disease) and other column represent markers.
  y <- as.numeric(as.matrix(table[,1]))
  X <- as.matrix(table[,-1])
  ## checking arguments
  Xy_perm <- check_args(y, X, perm=100)
  y <- Xy_perm$y
  X <- Xy_perm$X
  
  ## Number of individuals:
  #N <- nrow(X)
  ## Minor allele frequency (MAF):
  #MAF <- colMeans(X, na.rm=T)/2
  ## Number of variants with frequency < maf:
  #rare.maf <- MAF < maf
  #rare <- sum(rare.maf)
  #
  ## Collapsing procedure:
  #if (rare <= 1) {
  #  X.new <- X # We do not need to collapse
  #} else {
  #  X.collaps <- rowSums(X[, rare.maf], na.rm=T)
  #  X.collaps[X.collaps != 0] <- 1
  #  X.new <- cbind(X[, !rare.maf], X.collaps)
  #}
  
  cmc <- cmc_method(y, X, method)
  cmc.pval <- cmc$p.value
  #---------------------------#
  
  ## results
  name <- "CMC Test"
  arg.spec <- c(sum(y), length(y)-sum(y), ncol(X))
  arg.spec <- as.character(arg.spec)
  names(arg.spec) <- c("cases", "controls", "variants")
  res <- list(cmc.pval = cmc.pval,
              args = arg.spec, 
              name = name)
  return(res)
}