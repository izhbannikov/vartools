\name{vartools-package}
\alias{vartools-package}
\alias{vartools}
\docType{package}
\title{
Variant Association tools
}
\description{
\bold{vartools} is an R-package for the analysis of rare variants and currently has the following methods: CMC, KBAC,  RVE, VT, RareCover, RBT, C-alpha, etc.

\bold{Installation} \cr git clone https://github.com/izhbannikov/vartools.git \cr R CMD INSTALL vartools
}


\details{
\tabular{ll}{
Package: \tab vartools\cr
Type: \tab Package\cr
Version: \tab 1.3.0\cr
Date: \tab 2015-08-21\cr
License: \tab GPL (>= 2) \cr
Project URL: \tab \url{https://github.com/izhbannikov/vartools}
}

}
\author{
Ilya Y. Zhbannikov

Maintainer: Ilya Y. Zhbannikov | i.zhbannikov@mail.ru

}
\references{
}

\keyword{ package }


\examples{
#CMC test
# Load the package
library(vartools)
?cmc

pgdata <- as.matrix(read.table(system.file("extdata","phengen.dat",package="vartools"), as.is=T, skip = 1))
cmc.test <- cmc(table=pgdata)
print(cmc.test)

#KBAC test
?kbac

alpha <- 0.05
num.permutation <- 3000
quiet <- 1
alternative <- 1
maf.upper.bound <- 0.05
casectrl.dat <- read.table(system.file("extdata","phengen.dat",package="vartools"), skip = 1) 
kbac.test <- kbac(table=casectrl.dat, alpha, num.permutation, quiet, maf.upper.bound, alternative)
print(kbac.test)

}