 \name{rve}
 \alias{rve}
 \title{RVE statistic implementation}
 \description{
 This program implements the RVE statistic in [Cohen JC, Kiss RS, Pertsemlidis A, Marcel YL, McPherson R, et al. (2004)]. The RVE was first introduced in the analysis of seuence data from Dallas Heart Study. Performs a Kruskal-Wallis rank sum test.


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

## Hollander & Wolfe (1973), 116.
## Mucociliary efficiency from the rate of removal of dust in normal
##  subjects, subjects with obstructive airway disease, and subjects
##  with asbestosis.
x <- c(2.9, 3.0, 2.5, 2.6, 3.2) # normal subjects
y <- c(3.8, 2.7, 4.0, 2.4)      # with obstructive airway disease
z <- c(2.8, 3.4, 3.7, 2.2, 2.0) # with asbestosis
kruskal.test(list(x, y, z))
## Equivalently,
x <- c(x, y, z)
g <- factor(rep(1:3, c(5, 4, 5)),
            labels = c("Normal subjects",
                       "Subjects with obstructive airway disease",
                       "Subjects with asbestosis"))
kruskal.test(x, g)

## Formula interface.
require(graphics)
boxplot(Ozone ~ Month, data = airquality)
kruskal.test(Ozone ~ Month, data = airquality)


}
