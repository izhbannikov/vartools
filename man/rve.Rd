 \name{rve}
 \alias{rve}
 \title{RVE statistic implementation}
 \description{
 This program implements the RVE statistic in [Cohen JC, Kiss RS, Pertsemlidis A, Marcel YL, McPherson R, et al. (2004)]. The RVE was first introduced in the analysis of seuence data from Dallas Heart Study. Performs a one-sided Wilcoxon's two-sample test.


}
\usage{
rve(x,y)
}

\arguments{
      \item{x}{numeric vector of data values, DO NOT allow for missing data (will be omitted). }
      \item{y}{an optional numeric vector of data values, missing data will be omitted.}
}

\value{
	A list with class "htest" containing the following components:
\item{statistic}{the value of the test statistic with a name describing it.}
\item{parameter}{the parameter(s) for the exact distribution of the test statistic.}
\item{p.value}{the p-value for the test.}
\item{null.value}{the location parameter mu.}
\item{alternative}{a character string describing the alternative hypothesis.}
\item{method}{the type of test applied.}
\item{data.name}{a character string giving the names of the data.}
\item{conf.int}{a confidence interval for the location parameter. (Only present if argument conf.int = TRUE.)}
\item{estimate}{an estimate of the location parameter. (Only present if argument conf.int = TRUE.)}


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

x <- c(0.80, 0.83, 1.89, 1.04, 1.45, 1.38, 1.91, 1.64, 0.73, 1.46)
y <- c(1.15, 0.88, 0.90, 0.74, 1.21)
rve.stat <- rve(x, y)
rve.stat

}
