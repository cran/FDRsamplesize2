% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/samplesize.R
\name{find.sample.size}
\alias{find.sample.size}
\title{Determines the sample size needed to achieve the desired FDR and average power}
\usage{
find.sample.size(alpha, pwr, avepow.func, n0 = 3, n1 = 6, max.its = 50, ...)
}
\arguments{
\item{alpha}{the fixed p-value threshold  (scalar numeric)}

\item{pwr}{desired average power (scalar numeric)}

\item{avepow.func}{an R function to compute average power}

\item{n0}{lower limit for initial sample size range}

\item{n1}{upper limit for initial sample size range}

\item{max.its}{maximum number of iterations}

\item{...}{additional arguments to average power function}
}
\value{
A list with the following components:
\item{n}{a sample size estimate}
\item{computed.avepow}{average power}
\item{desired.avepow}{desired average power}
\item{alpha}{fixed p-value threshold for multiple testing procedure}
\item{n.its}{number of iteration}
\item{max.its}{maximum number of iteration, default is 50}
\item{n0}{lower limit for initial sample size range}
\item{n1}{upper limit for initial sample size range}
}
\description{
Determines the sample size needed to achieve the desired FDR and average power by given the proportion of true null hypothesis.
}
\note{
For the test with power calculation based on asymptotic normal approximation, we suggest checking \code{FDRsamplesize2} calculation by simulation.
}
\examples{
#Here, calculating the sample size for the study involving many sign tests
average.power.signtest;
p.adj = 0.001;
p = rep(c(0.8,0.5), c(100,9900));
find.sample.size(alpha = p.adj, pwr = 0.8, avepow.func = average.power.signtest, p = p)
}
