 \name{calpha}
 \alias{calpha}
 \title{C(alpha) statistical test}
 \description{
 C(alpha) (Neale at al., 2011) is a statistical test, developed for disease traits, in order to test for the hypothesis of rare variants disease association. This test has an assumption that rare variants, which were observed in cases and controls and are a mix of various types of variants: phenotypically deleterious, protective and neutral variants.
}
\usage{
calpha(table, permutations=NULL)
}

\arguments{
      \item{table}{a numeric matrix with first column having disease status '0' or '1' and the rest columns codes the locus genotype as '0', '1', and '2'.}
      \item{permutations}{positive integer that defines the number of permutations. In permutation test, the distribution of the test statistic under the null hypothesis is obtained by calculating all possible values of the test statistic under rearrangements of the labels on the observed data points.}
}

\value{
A list with the following elements:
\item{calpha.stat}{C(alpha) statistic}
\item{asym.pval}{asymptotic p-value, distributed as Chi-square, with parameter degrees of freedom \code{df=1}}
\item{perm.pval}{permuted p-value; only when \code{perm} is used}
\item{args}{descriptive information with number of controls, cases,variants, and permutations}
\item{name}{name of the statistic}
}

\details{
    ....
}

\author{}

\references{
Benjamin M. Neale, Manuel A. Rivas, Benjamin F. Voight, David Altshuler, Bernie Devlin, Marju Orho-Melander, Sekar Kathiresan, Shaun M. Purcell, Kathryn Roeder and Mark J. Daly (2011), Testing for an Unusual Distribution of Rare Variants. PLoS Genetics doi:10.1371/journal.pgen.1001322. http://dx.plos.org/10.1371/journal.pgen.1001322
}

\examples{
# Load the package

library(vartools)
?calpha

casectrl.dat <- read.table(system.file("extdata","phengen.dat",package="vartools"), skip = 1)    
calpha.stat <- calpha(casectrl.dat)
calpha.stat
}
