 \name{cmc}
 \alias{cmc}
 \title{CMC statistic implementation}
 \description{
 This program implements the CMC statistic in [Liu and Leal 2008]. It carries out case-control association testing for rare variants for whole exome association studies. Briefly, consider a gene of length n which harbors m rare variants. Genotype on the m variant sites & the disease status (case/control) are known for each individual. The program takes as input the m-site genotype and disease status (case/control) data files, and computes a p-value indicating the significance of association. In order to speed up permutation testing we use an "adaptive" approach to obtain p-values.
}

\usage{
cmc(table, method = "fisher")
}

\arguments{
      \item{table}{a numeric matrix with first column having disease status '0' or '1' and the rest columns codes the locus genotype as '0', '1', and '2'. DO NOT allow for missing data. }
      \item{method}{statistical method, Fisher test used by default}
}

\value{
	\item{pvalue}{the p-value of test.}
}

\details{
    ....
}

\author{Ilya Y. Zhbannikov | i.zhbannikov@mail.ru}

\references{
Liu DJ,  Leal SM, 2008 Methods for Detecting Associations
with Rare Variants for Common Diseases: Application to Analysis of Sequence Data. The American Journal of Human Genetics, 83, DOI10.1016/j.ajhg.2008.06.024  
}

\examples{
# Load the package
library(vartools)
?cmc

pgdata <- as.matrix(read.table(system.file("extdata","phengen.dat",package="vartools"), as.is=T, skip = 1))
cmc.pvalue <- cmc(table=pgdata)
print(cmc.pvalue)

}
