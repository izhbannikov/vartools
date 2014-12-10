 \name{rarecover}
 \alias{rarecover}
 \title{RareCover statistical test}
 \description{
 Testing for rare variants with RareCover algorithm. This algorithm is similar to CMC, meaning that it follows its collapsing strategy, but uses greedy algorithm to find an optimized combination of variants in a loci for which its association signal is strongest.
}
\usage{
rarecover(table, maf=0.05, dif=0.5, perm=250)
}

\arguments{
      \item{table}{a numeric matrix with first column having disease status '0' or '1' and the rest columns codes the locus genotype as '0', '1', and '2'.}
      \item{maf}{numeric value indicating the minor allele frequency threshold for rare variants }
      \item{dif}{numeric value between 0 and 1 as a threshold for the decision criterion in the RareCover algorithm}
      \item{perm}{positive integer that defines the number of permutations}
}

\value{
A list with the following elements:
\item{rc.stat}{RareCover statistic}
\item{perm.pval} {permuted p-value}
\item{set}{set of selected variants}
\item{args}{descriptive information with number of controls, cases,variants, and permutations}
\item{name} {name of the statistic}
}

\details{
    ....
}

\author{}

\references{
Bhatia G, Bansal V, Harismendy O, Schork NJ, Topol EJ, Frazer K, Bafna V (2010) A Covering Method for Detecting Genetic Associations between Rare Variants and Common Phenotypes. \emph{PLoS Computational Biology}, \bold{6(10)}: e1000954
}

\examples{
# Load the package

library(vartools)
?rarecover

casectrl.dat <- read.table(system.file("extdata","phengen.dat",package="vartools"), skip = 1)
rarecover.stat <- rarecover(casectrl.dat)
rarecover.stat
}
