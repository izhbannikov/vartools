cmc <- function(filename, method="fisher") {
  pgdata <- as.matrix(read.table(filename, as.is=T, skip = 1))
  y <- pgdata[,1];
  x <- ((rowSums(as.matrix(pgdata[,-1])))>0);
  # Contigency table:
  m <- matrix(nrow=2,ncol=2);
  m[1,1] <- sum(x==1 & y==1);
  m[1,2] <- sum(x==0 & y==1);
  m[2,1] <- sum(x==1 & y==0);
  m[2,2] <- sum(x==0 & y==0);
  print(m);
  if (method == "fisher") {
    stat <- fisher.test(m)$p.value; # Fisher's exact test
  }
  if (method == "chisq") {
    stat <- chisq.test(m)$p.value;
  }
  
  return(stat);
}