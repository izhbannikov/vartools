cmc <- function(table, method="fisher") {
  # Input table consists of columns and rows. Rows are individuals and columns are:
  # First column represents a status (disease / no disease) and other column represent markers.
  pgdata <- as.matrix(read.table(table, as.is=T, skip = 1))
  markers <- pgdata[,1]; # Markers
  status <- ((rowSums(as.matrix(pgdata[,-1])))>0);
  
  # Contingency table:
  ct <- matrix(nrow=2,ncol=2);
  ct[1,1] <- sum(markers==1 & status==1);
  ct[1,2] <- sum(markers==0 & status==1);
  ct[2,1] <- sum(markers==1 & status==0);
  ct[2,2] <- sum(markers==0 & status==0);
  
  print(ct); # Let's take a look at the data
  
  if (method == "fisher") {
    stat <- fisher.test(ct)$p.value; # Fisher's exact test
  }
  if (method == "chisq") {
    stat <- chisq.test(ct)$p.value; # Should this test work?
  }
  
  stat;
}