 \name{rve}
 \alias{rve}
 \title{RVE statistic implementation}
 \description{
 This program implements the RVE statistic in [Cohen JC, Kiss RS, Pertsemlidis A, Marcel YL, McPherson R, et al. (2004)]. The RVE was first introduced in the analysis of seuence data from Dallas Heart Study. It compares frequency of carriers of rare variants that are found exclusively in cases or controls using Fisher's exact test.
}
\usage{
rve(table)
}

\arguments{
      \item{table}{a contingency table, DO NOT allow for missing data. }
}

\value{
	\item{pvalue}{the p-value of test.}
}

\details{
    ....
}

\author{Ilya Y. Zhbannikov | i.zhbannikov@mail.ru}

\references{
Cohen JC, Kiss RS, Pertsemlidis A, Marcel YL, McPherson R, et al. (2004) Multiple rare alleles contribute to low plasma levels of HDL cholesterol. Science 305: 869-872.
}

\examples{
# Load the package

library(vartools)
?rve

c.table <- matrix(c(13, 115, 3, 125),
                  nrow = 2,
                  dimnames = list(c("NS", "NO_NS"),
                                  c("L", "H")))
                                  
rve.pvalue <- rve(c.table)
print(rve.pvalue)

}
