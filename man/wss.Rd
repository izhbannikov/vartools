 \name{wss}
 \alias{wss}
 \title{WSS statistical test}
 \description{
 WSS (Madsen and Browning, 2009) is a statistical test, developed for disease traits, in order to test for the hypothesis of rare variants disease association. WSS test introduces the method of assigning "weights" to rare variants found in a genetic region before they are collapsed.
}
\usage{
wss(table, perm=100)
}

\arguments{
      \item{table}{a numeric matrix with first column having disease status '0' or '1' and the rest columns codes the locus genotype as '0', '1', and '2'.}
      \item{perm}{positive integer that defines the number of permutations. In permutation test, the distribution of the test statistic under the null hypothesis is obtained by calculating all possible values of the test statistic under rearrangements of the labels on the observed data points.}
}

\value{
A list with the following elements:
\item{wss.stat}{WSS statistic}
\item{asym.pval}{asymptotic p-value}
\item{perm.pval}{permuted p-value; only when \code{permutations} is used}
\item{args}{descriptive information with number of controls, cases,variants, and permutations}
\item{name}{name of the statistic}
}

\details{
    ....
}

\author{}

\references{
Bo Eskerod Madsen and Sharon R. Browning (2009) A Groupwise Association Test for Rare Mutations Using a Weighted Sum Statistic. PLoS Genetics doi:10.1371/journal.pgen.1000384. http://dx.plos.org/10.1371/journal.pgen.1000384
}

\examples{
# Load the package

library(vartools)
?wss

casectrl.dat <- read.table(system.file("extdata","phengen.dat",package="vartools"), skip = 1)    
wss.stat <- wss(casectrl.dat)
wss.stat
}
