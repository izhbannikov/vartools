\name{vartools-package}
\alias{vartools-package}
\alias{vartools}
\docType{package}
\title{
Variant Association tools
}
\description{
\bold{vartools} is an R-package for the analysis of variants and currently has two methods: \bold{CMC} and \bold{KBAC}.
}
\details{
\tabular{ll}{
Package: \tab vartools\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2014-11-27\cr
License: \tab GPL (>= 2) \cr
}

Package contains the following functions:

\bold{CMC} - implements the CMC statistic in [Liu and Leal 2008]. It carries out case-control association testing for rare variants for whole exome association studies.

\bold{KBAC} - implements the KBAC statistic in [Liu and Leal 2010].

}
\author{
Ilya Y. Zhbannikov

Maintainer: Ilya Y. Zhbannikov | i.zhbannikov@mail.ru

}
\references{
Liu DJ,  Leal SM, 2010 A Novel Adaptive Method for the Analysis of Next-Generation Sequencing Data to Detect Complex Trait Associations with Rare Variants Due to Gene Main Effects and Interactions. PLoS Genet 6(10): e1001156. doi:10.1371/journal.pgen.1001156 

Liu DJ,  Leal SM, 2008 Methods for Detecting Associations
with Rare Variants for Common Diseases: Application to Analysis of Sequence Data. The American Journal of Human Genetics, 83, DOI10.1016/j.ajhg.2008.06.024 
}

\keyword{ package }

\examples{
#CMC test
# Load the package
library(vartools)
?cmc

cmc.pvalue <- cmc(table=system.file("extdata","phengen.dat",package="vartools"))
print(cmc.pvalue)
summary(cmc.pvalue)

#KBAC test
alpha <- 0.05
num.permutation <- 3000
quiet <- 1
alternative <- 1
maf.upper.bound <- 0.05
kbac.pvalue <- kbac(table=system.file("extdata","phengen.dat",package="vartools"), alpha, num.permutation, quiet, maf.upper.bound, alternative)
print(kbac.pvalue)
summary(kbac.pvalue)
}