 \name{vt}
 \alias{vt}
 \title{Variable Threshold statistical test}
 \description{
 This is implementation for the Variable Threshold statistical test. It assigns different weights to rare variants base on their frequencies.
}
\usage{
vt(table, maf=0.05, perm=50)
}

\arguments{
      \item{table}{a numeric matrix with first column having disease status '0' or '1' and the rest columns codes the locus genotype as '0', '1', and '2'.}
      \item{maf}{numeric value indicating the minor allele frequency threshold for rare variants (must be a positive number between 0 and 1, \code{maf=0.05} by default).}
      \item{perm}{positive integer that defines the number of permutations, \code{perm=50} by default. In permutation test, the distribution of the test statistic under the null hypothesis is obtained by calculating all possible values of the test statistic under rearrangements of the labels on the observed data points.}
}

\value{
A list with the following elements:
\item{vt.stat}{VT statistic}
\item{perm.pval}{permuted p-value}
\item{args}{descriptive information with number of controls, cases,variants, and permutations}
\item{name}{name of the statistic}
}

\details{
    ....
}

\author{}

\references{
Price AL, Kryukov GV, de Bakker PIW, Purcell SM, Staples J, Wei LJ, Sunyaev SR (2010) Pooled Association Tests for Rare Variants in Exon-Sequencing Studies. \emph{The American Journal of Human Genetics}, \bold{86}: 832-838
}

\examples{
# Load the package
library(vartools)
?vt

casectrl.dat <- read.table(system.file("extdata","phengen2.dat",package="vartools"), skip = 1)
vt.stat <- vt(casectrl.dat)
vt.stat

}
