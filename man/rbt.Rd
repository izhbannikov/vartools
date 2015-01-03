 \name{rbt}
 \alias{rbt}
 \title{Replication-Based statistical test}
 \description{
  Implementation of the replication base test (RBT) proposed by Ionita-Laza et al., 2011. 
Replication means that it computes evidences to reject each of the two hypothesis:
(1) Deleterious rare variants are enriched in cases; (2) Protective rare variants are enriched in controls.
}
\usage{
rbt(table, perm=150)
}

\arguments{
      \item{table}{a numeric matrix with first column having disease status '0' or '1' and the rest columns codes the locus genotype as '0', '1', and '2'.}
      \item{perm}{positive integer that defines the number of permutations, \code{150} by default. In permutation test, the distribution of the test statistic under the null hypothesis is obtained by calculating all possible values of the test statistic under rearrangements of the labels on the observed data points.}
}

\value{
A list with the following elements:
\item{rbt.stat}{RBT statistic}
\item{perm.pval}{permuted p-value}
\item{args}{descriptive information with number of controls, cases,variants, and permutations}
\item{name}{name of the statistic}
}

\details{
    ....
}

\author{}

\references{
Ionita-Laza I, Buxbaum JD, Laird NM, Lange C (2011) A New Testing Strategy to Identify Rare Variants with Either risk or Protective Effects on Disease. \emph{PLoS Genetics}, \bold{7(2)}: e1001289
}

\examples{
# Load the package

library(vartools)
?rbt

casectrl.dat <- read.table(system.file("extdata","phengen.dat",package="vartools"), skip = 1)
rbt.stat <- rbt(casectrl.dat)
rbt.stat
}
