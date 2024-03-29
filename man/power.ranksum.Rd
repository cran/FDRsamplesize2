% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/power.R
\name{power.ranksum}
\alias{power.ranksum}
\title{Compute power of the rank-sum test}
\usage{
power.ranksum(n, alpha, p)
}
\arguments{
\item{n}{sample size (scalar)}

\item{alpha}{p-value threshold (scalar)}

\item{p}{Pr(Y>X), as in Noether (JASA 1987)}
}
\value{
Vector of power estimates for two-sided tests
}
\description{
Compute power of rank-sum test; Uses formula of Noether (JASA 1987), which is based on asymptotic normal approximation.
}
\details{
In most applications, the null effect size will be designated by p = 0.5
}
\examples{
p = rep(c(0.8,0.5),c(100,900))
res = power.ranksum(n = 50, alpha = 0.5, p=p)
}
\references{
Noether, Gottfried E (1987) Sample size determination for some common nonparametric tests. Journal of the American Statistical Association, 82:645-647.
}
