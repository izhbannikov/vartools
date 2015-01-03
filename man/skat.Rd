 \name{skat}
 \alias{skat}
 \title{SKAT statistical test}
 \description{
 SKAT (Wu at al., 2010, 2011) is a statistical test, developed for disease traits, in order to test for the hypothesis of rare variants disease association. This is a score-based variance-component test.
}
\usage{
skat(table, kernel="linear", weights=NULL, a=1, b=25, permutations=NULL)
}

\arguments{
      \item{table}{a numeric matrix with first column having disease status '0' or '1' and the rest columns codes the locus genotype as '0', '1', and '2'.}
      \item{kernel}{Kernel function used. By default, a linear function is used. Other available kernels are: \code{"wlinear"}, \code{"quadratic"}, \code{"IBS"}, \code{"wIBS"}, \code{"twowayx"}} 
      \item{weights}{a set of numeric weights for genetic variants.}
      \item{a}{a positive numeric value for the parameter \code{a} in the Beta distribution (\code{a=1} by default)}
      \item{b}{a positive numeric value for the parameter \code{b} in the Beta distribution (\code{b=1} by default)}
      \item{permutations}{positive integer that defines the number of permutations. In permutation test, the distribution of the test statistic under the null hypothesis is obtained by calculating all possible values of the test statistic under rearrangements of the labels on the observed data points.}
}

\value{
A list with the following elements:
\item{skat.stat}{SKAT statistic}
\item{asym.pval}{asymptotic p-value, distributed as Chi-square, with parameter degrees of freedom \code{df=1}}
\item{perm.pval}{permuted p-value; only when \code{permutations} is used}
\item{args}{descriptive information with number of controls, cases,variants, and permutations}
\item{name}{name of the statistic}
}

\details{
    ....
}

\author{}

\references{
Wu MC, Kraft P, Epstein MP, Taylor DM, Chanock SJ, Hunter DJ, Lin X (2010) Powerful SNP-Set Analysis for Case-Control Genome-wide Association Studies. \emph{The American Journal of Human Genetics}, \bold{86}: 929-942
\cr
Wu MC, Lee S, Cai T, Li Y, Boehnke M, Lin X (2011) Rare-Variant Association Testing for Sequencing Data with the Sequence Kernel Association Test. \emph{The American Journal of Human Genetics}, \bold{89}: 82-93
}

\examples{
# Load the package

library(vartools)
?skat

casectrl.dat <- read.table(system.file("extdata","phengen.dat",package="vartools"), skip = 1)    
skat.stat <- skat(casectrl.dat)
skat.stat
}
