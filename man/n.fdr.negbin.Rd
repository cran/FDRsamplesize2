% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/samplesize.R
\name{n.fdr.negbin}
\alias{n.fdr.negbin}
\title{Sample size calculation for Negative Binomial data}
\usage{
n.fdr.negbin(fdr, pwr, log.fc, mu, sig, pi0.hat = "BH")
}
\arguments{
\item{fdr}{desired FDR (scalar numeric)}

\item{pwr}{desired average power (scalar numeric)}

\item{log.fc}{log fold-change (vector), usual null hypothesis is log.fc=0}

\item{mu}{read depth per gene (vector, same length as log.fc)}

\item{sig}{coefficient of variation (CV) per gene (vector, same length as log.fc)}

\item{pi0.hat}{method to estimate proportion \code{pi0} of tests with true null, including: "HH" (p-value histogram height), "HM" (p-value histogram mean), "BH" (Benjamini & Hochberg 1995), "Jung" (Jung 2005)}
}
\value{
A list with the following components:
\item{n}{per-group sample size estimate}
\item{computed.avepow}{average power}
\item{desired.avepow}{desired average power}
\item{desired.fdr}{desired FDR}
\item{input.pi0}{proportion of tests with a true null hypothesis}
\item{alpha}{fixed p-value threshold for multiple testing procedure}
\item{n.its}{number of iteration}
\item{max.its}{maximum number of iteration, default is 50}
\item{n0}{lower limit for initial sample size range}
\item{n1}{upper limit for initial sample size range}
}
\description{
Find the sample size needed to have a desired false discovery rate and average power for a large number of Negative Binomial comparisons.
}
\note{
For the test with power calculation based on asymptotic normal approximation, we suggest checking \code{FDRsamplesize2} calculation by simulation.
}
\examples{
logFC = log(rep(c(1,2),c(900,100)));
mu = rep(5,1000);
sig = rep(0.6,1000);
n.fdr.negbin(fdr = 0.1, pwr = 0.8, log.fc = logFC, mu = mu, sig = sig, pi0.hat = "BH")
}
\references{
SN Hart, TM Therneau, Y Zhang, GA Poland, and J-P Kocher (2013). Calculating Sample Size Estimates for RNA Sequencing Data. Journal of Computational Biology 20: 970-978.
}
