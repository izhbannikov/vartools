 \name{asum}
 \alias{asum}
 \title{aSum statistical test}
 \description{
 aSum, Adaptive Sum Test (Han and Pan, 2010) is a statistical test which utilizes the difference in direction of effects (protective or deleterious) of rare variants within the same genetic region analyzed by a rare variant association test.
}
\usage{
asum(table, perm=100)
}

\arguments{
      \item{table}{a numeric matrix with first column having disease status '0' or '1' and the rest columns codes the locus genotype as '0', '1', and '2'.}
      \item{perm}{positive integer that defines the number of permutations. In permutation test, the distribution of the test statistic under the null hypothesis is obtained by calculating all possible values of the test statistic under rearrangements of the labels on the observed data points.}
}

\value{
A list with the following elements:
\item{asum.stat}{aSum statistic}
\item{perm.pval}{permuted p-value; only when \code{permutations} is used}
\item{args}{descriptive information with number of controls, cases,variants, and permutations}
\item{name}{name of the statistic}
}

\details{
    ....
}

\author{}

\references{
Fang Han and Wei Pan (2010) A Data-Adaptive Sum Test for Disease Association with Multiple Common or Rare Variants. Human Heredity doi:10.1159/000288704. http://www.karger.com/doi/10.1159/000288704
}

\examples{
# Load the package

library(vartools)
?asum

casectrl.dat <- read.table(system.file("extdata","phengen.dat",package="vartools"), skip = 1)    
asum.stat <- asum(casectrl.dat)
asum.stat
}
