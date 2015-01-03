# These are Kernels (functions) used in association tests
kernel_IBS <- function(Z, n, p) {
    K  <-  diag(1, n, n)	
    aux  <-  .C("kernel_IBS", as.integer(as.vector(t(Z))), as.integer(n), as.integer(p), as.double(as.vector(K)))[[4]]
    matrix(aux, nrow=n)
}

kernel_wIBS <- function(Z, n, p, weights) {
    K <- diag(1, n, n)
    aux <- .C("kernel_wIBS", as.integer(as.vector(t(Z))), as.integer(n), as.integer(p), 
             as.double(weights), as.double(as.vector(K)))[[5]]
    matrix(aux, nrow=n)
}

kernel_twowayx <- function(Z, n, p) {
    K <- matrix(0, n, n)	
    aux <- .C("kernel_twowayx", as.integer(as.vector(t(Z))), as.integer(n), as.integer(p), as.double(as.vector(K)))[[4]]
    matrix(aux, nrow=n)
}


