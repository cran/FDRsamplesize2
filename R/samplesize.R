#' @title Sample size calculation for t-tests for non-zero correlation
#' @description
#' Find the sample size needed to have a desired false discovery rate and average power for a large number of t-tests for non-zero correlation.
#' @param fdr desired FDR (scalar numeric)
#' @param pwr desired average power (scalar numeric)
#' @param rho  population correlation coefficient (vector)
#' @param pi0.hat method to estimate proportion \code{pi0} of tests with true null, including: "HH" (p-value histogram height), "HM" (p-value histogram mean), "BH" (Benjamini & Hochberg 1995), "Jung" (Jung 2005)
#' @return  A list with the following components:
#' \item{n}{sample size estimate}
#' \item{computed.avepow}{average power}
#' \item{desired.avepow}{desired average power}
#' \item{desired.fdr}{desired FDR}
#' \item{input.pi0}{proportion of tests with a true null hypothesis}
#' \item{alpha}{fixed p-value threshold for multiple testing procedure}
#' \item{n.its}{number of iteration}
#' \item{max.its}{maximum number of iteration, default is 50}
#' \item{n0}{lower limit for initial sample size range}
#' \item{n1}{upper limit for initial sample size range}
#' @examples
#' rho = rep(c(0.3,0),c(100,900));
#' n.fdr.tcorr(fdr = 0.1, pwr = 0.8, rho = rho, pi0.hat="BH")
#' @export
n.fdr.tcorr=function(fdr,pwr,rho,pi0.hat="BH")
{
  pi0=mean(rho==0)
  a=alpha.power.fdr(fdr,pwr,pi0,pi0.hat)
  res=find.sample.size(alpha=a,pwr=pwr,avepow.func=average.power.tcorr, rho=rho)
  res$desired.fdr=fdr
  res$input.pi0=pi0
  res.names=c("n",
              "computed.avepow",
              "desired.avepow",
              "desired.fdr",
              "input.pi0",
              "alpha",
              "n0","n1",
              "n.its",
              "max.its")
  res=res[res.names]
  return(res)

}



#' @title Sample size calculation for comparing two proportions
#' @description
#' Find the sample size needed to have a desired false discovery rate and average power for a large number of two-group comparisons using the two proportion z-test.
#' @param fdr desired FDR (scalar numeric)
#' @param pwr desired average power (scalar numeric)
#' @param p1  probability in one group (vector)
#' @param p2  probability in other group (vector)
#' @param alternative one- or two-sided test
#' @param pi0.hat method to estimate proportion \code{pi0} of tests with true null, including: "HH" (p-value histogram height), "HM" (p-value histogram mean), "BH" (Benjamini & Hochberg 1995), "Jung" (Jung 2005)
#' @return  A list with the following components:
#' \item{n}{per-group sample size estimate}
#' \item{computed.avepow}{average power}
#' \item{desired.avepow}{desired average power}
#' \item{desired.fdr}{desired FDR}
#' \item{input.pi0}{proportion of tests with a true null hypothesis}
#' \item{alpha}{fixed p-value threshold for multiple testing procedure}
#' \item{n.its}{number of iteration}
#' \item{max.its}{maximum number of iteration, default is 50}
#' \item{n0}{lower limit for initial sample size range}
#' \item{n1}{upper limit for initial sample size range}
#' @note For the test with power calculation based on asymptotic normal approximation, we suggest checking \code{FDRsamplesize2} calculation by simulation.
#' @examples
#' set.seed(1234);
#' p1 = sample(seq(0,0.5,0.1),40,replace = TRUE);
#' p2 = sample(seq(0.5,1,0.1),40,replace = TRUE);
#' n.fdr.twoprop(fdr = 0.1, pwr = 0.8, p1 = p1, p2 = p2, alternative = "two.sided", pi0.hat = "BH")
#' @export
n.fdr.twoprop <- function(fdr, pwr, p1, p2, alternative = "two.sided", pi0.hat = "BH") {
  pi0 <- mean(p1 == p2)
  a <- alpha.power.fdr(fdr, pwr, pi0, pi0.hat)
  res <- find.sample.size(alpha = a, pwr = pwr, avepow.func = average.power.twoprop, p1 = p1, p2 = p2, alternative = alternative)
  res$desired.fdr=fdr
  res$input.pi0=pi0
  res.names=c("n",
              "computed.avepow",
              "desired.avepow",
              "desired.fdr",
              "input.pi0",
              "alpha",
              "n0","n1",
              "n.its",
              "max.its")
  res=res[res.names]
  return(res)
}


#' @title Sample size calculation for Fisher's Exact tests
#' @description
#' Find the sample size needed to have a desired false discovery rate and average power for a large number of Fisher's exact tests.
#' @param fdr desired FDR (scalar numeric)
#' @param pwr desired average power (scalar numeric)
#' @param p1  probability in one group (vector)
#' @param p2  probability in other group (vector)
#' @param alternative one- or two-sided test
#' @param pi0.hat method to estimate proportion \code{pi0} of tests with true null, including: "HH" (p-value histogram height) , "HM" (p-value histogram mean), "BH" (Benjamini & Hochberg 1995), "Jung" (Jung 2005)
#' @return  A list with the following components:
#' \item{n}{per-group sample size estimate}
#' \item{computed.avepow}{average power}
#' \item{desired.avepow}{desired average power}
#' \item{desired.fdr}{desired FDR}
#' \item{input.pi0}{proportion of tests with a true null hypothesis}
#' \item{alpha}{fixed p-value threshold for multiple testing procedure}
#' \item{n.its}{number of iteration}
#' \item{max.its}{maximum number of iteration, default is 50}
#' \item{n0}{lower limit for initial sample size range}
#' \item{n1}{upper limit for initial sample size range}
#' @examples
#' set.seed(1234);
#' p1 = sample(seq(0,0.5,0.1),10,replace = TRUE);
#' p2 = sample(seq(0.5,1,0.1),10,replace = TRUE);
#' n.fdr.fisher(fdr = 0.1, pwr = 0.8, p1 = p1, p2 = p2, alternative = "two.sided", pi0.hat = "BH")
#' @export
n.fdr.fisher <- function(fdr, pwr, p1, p2, alternative = "two.sided", pi0.hat = "BH") {
  pi0 <- mean(p1 == p2)
  a <- alpha.power.fdr(fdr, pwr, pi0, pi0.hat)
  res <- find.sample.size(alpha = a, pwr = pwr, avepow.func = average.power.fisher, p1 = p1, p2 = p2, alternative = alternative)
  res$desired.fdr=fdr
  res$input.pi0=pi0
  res.names=c("n",
              "computed.avepow",
              "desired.avepow",
              "desired.fdr",
              "input.pi0",
              "alpha",
              "n0","n1",
              "n.its",
              "max.its")
  res=res[res.names]
  return(res)
}


#' @title Sample size calculation for Poisson data
#' @description
#' Find the sample size needed to have a desired false discovery rate and average power for a large number of two-group comparisons under Poisson distribution.
#' @param fdr desired FDR (scalar numeric)
#' @param pwr desired average power (scalar numeric)
#' @param rho fold-change, usual null hypothesis is that rho=1 (vector)
#' @param mu0 average count in control group (vector)
#' @param w ratio of the total number of reads mapped between the two groups
#' @param type type of test: "w" for Wald, "s" for score, "lw" for log-transformed Wald, "ls" for log-transformed score.
#' @param pi0.hat method to estimate proportion \code{pi0} of tests with true null, including: "HH" (p-value histogram height) , "HM" (p-value histogram mean), "BH" (Benjamini & Hochberg 1995), "Jung" (Jung 2005)
#' @return  A list with the following components:
#' \item{n}{per-group sample size estimate}
#' \item{computed.avepow}{average power}
#' \item{desired.avepow}{desired average power}
#' \item{desired.fdr}{desired FDR}
#' \item{input.pi0}{proportion of tests with a true null hypothesis}
#' \item{alpha}{fixed p-value threshold for multiple testing procedure}
#' \item{n.its}{number of iteration}
#' \item{max.its}{maximum number of iteration, default is 50}
#' \item{n0}{lower limit for initial sample size range}
#' \item{n1}{upper limit for initial sample size range}
#' @references C-I Li, P-F Su, Y Guo, and Y Shyr (2013). Sample size calculation for differential expression analysis of RNA-seq data under Poisson distribution. Int J Comput Biol Drug Des 6(4).<doi:10.1504/IJCBDD.2013.056830>
#' @examples
#' rho = rep(c(1,1.25),c(900,100));
#' mu0 = rep(5,1000);
#' w = rep(0.5,1000);
#' n.fdr.poisson(fdr = 0.1, pwr = 0.8, rho = rho, mu0 = mu0, w = w, type = "w", pi0.hat = "BH")
#' @export
n.fdr.poisson <- function(fdr, pwr, rho, mu0, w, type, pi0.hat = "BH") {
  pi0 <- mean(rho == 1)
  a <- alpha.power.fdr(fdr, pwr, pi0, pi0.hat)
  res <- find.sample.size(alpha = a, pwr = pwr, avepow.func = average.power.li, rho = rho, mu0 = mu0, w = w, type = type)
  res$desired.fdr=fdr
  res$input.pi0=pi0
  res.names=c("n",
              "computed.avepow",
              "desired.avepow",
              "desired.fdr",
              "input.pi0",
              "alpha",
              "n0","n1",
              "n.its",
              "max.its")
  res=res[res.names]
  return(res)
}


#' @title Sample size calculation for signed-rank tests
#' @description
#' Find the sample size needed to have a desired false discovery rate and average power for a large number of signed-rank tests.
#' @param fdr desired FDR (scalar numeric)
#' @param pwr desired average power (scalar numeric)
#' @param p1 Pr(X>0), as in Noether (JASA 1987)
#' @param p2 Pr(X+X'>0), as in Noether (JASA 1987)
#' @param pi0.hat method to estimate proportion \code{pi0} of tests with true null, including: "HH" (p-value histogram height) , "HM" (p-value histogram mean), "BH" (Benjamini & Hochberg 1995), "Jung" (Jung 2005)
#' @return  A list with the following components:
#' \item{n}{sample size estimate}
#' \item{computed.avepow}{average power}
#' \item{desired.avepow}{desired average power}
#' \item{desired.fdr}{desired FDR}
#' \item{input.pi0}{proportion of tests with a true null hypothesis}
#' \item{alpha}{fixed p-value threshold for multiple testing procedure}
#' \item{n.its}{number of iteration}
#' \item{max.its}{maximum number of iteration, default is 50}
#' \item{n0}{lower limit for initial sample size range}
#' \item{n1}{upper limit for initial sample size range}
#' @references Noether, Gottfried E (1987) Sample size determination for some common nonparametric tests. Journal of the American Statistical Association, 82:645-647.
#' @examples
#' p1 = rep(c(0.8,0.5),c(100,900));
#' p2 = rep(c(0.8,0.5),c(100,900));
#' n.fdr.signrank(fdr = 0.1, pwr = 0.8, p1 = p1, p2 = p2, pi0.hat = "BH")
#' @export
n.fdr.signrank <- function(fdr, pwr, p1, p2, pi0.hat = "BH") {
  pi0 <- length(intersect(which(p1 == 0.5),which(p2 == 0.5)))/length(p1)
  a <- alpha.power.fdr(fdr, pwr, pi0, pi0.hat)
  res <- find.sample.size(alpha = a, pwr = pwr, avepow.func = average.power.signrank, p1 = p1, p2 = p2)
  res$desired.fdr=fdr
  res$input.pi0=pi0
  res.names=c("n",
              "computed.avepow",
              "desired.avepow",
              "desired.fdr",
              "input.pi0",
              "alpha",
              "n0","n1",
              "n.its",
              "max.its")
  res=res[res.names]
  return(res)
}


#' @title Sample size calculation for sign tests
#' @description
#' Find the sample size needed to have a desired false discovery rate and average power for a large number of sign tests.
#' @param fdr desired FDR (scalar numeric)
#' @param pwr desired average power (scalar numeric)
#' @param p Pr(X>0), as in Noether (JASA 1987)
#' @param pi0.hat method to estimate proportion \code{pi0} of tests with true null, including: "HH" (p-value histogram height), "HM" (p-value histogram mean), "BH" (Benjamini & Hochberg 1995), "Jung" (Jung 2005)
#' @return  A list with the following components:
#' \item{n}{sample size estimate}
#' \item{computed.avepow}{average power}
#' \item{desired.avepow}{desired average power}
#' \item{desired.fdr}{desired FDR}
#' \item{input.pi0}{proportion of tests with a true null hypothesis}
#' \item{alpha}{fixed p-value threshold for multiple testing procedure}
#' \item{n.its}{number of iteration}
#' \item{max.its}{maximum number of iteration, default is 50}
#' \item{n0}{lower limit for initial sample size range}
#' \item{n1}{upper limit for initial sample size range}
#' @references Noether, Gottfried E (1987) Sample size determination for some common nonparametric tests. Journal of the American Statistical Association, 82:645-647.
#' @note For the test with power calculation based on asymptotic normal approximation, we suggest checking \code{FDRsamplesize2} calculation by simulation.
#' @examples
#' p = rep(c(0.8, 0.5), c(100, 900));
#' n.fdr.signtest(fdr = 0.1, pwr = 0.8, p = p, pi0.hat = "BH")
#' @export
n.fdr.signtest <- function(fdr, pwr, p, pi0.hat = "BH") {
  pi0 <- mean(p == 0.5)
  a <- alpha.power.fdr(fdr, pwr, pi0, pi0.hat)
  res <- find.sample.size(alpha = a, pwr = pwr, avepow.func = average.power.signtest, p = p)
  res$desired.fdr=fdr
  res$input.pi0=pi0
  res.names=c("n",
              "computed.avepow",
              "desired.avepow",
              "desired.fdr",
              "input.pi0",
              "alpha",
              "n0","n1",
              "n.its",
              "max.its")
  res=res[res.names]
  return(res)
}


#' @title Sample size calculation for rank-sum tests
#' @description
#' Find the sample size needed to have a desired false discovery rate and average power for a large number of rank-sum tests.
#' @param fdr desired FDR (scalar numeric)
#' @param pwr desired average power (scalar numeric)
#' @param p Pr(Y>X), as in Noether (JASA 1987)
#' @param pi0.hat method to estimate proportion \code{pi0} of tests with true null, including: "HH" (p-value histogram height), "HM" (p-value histogram mean), "BH" (Benjamini & Hochberg 1995), "Jung" (Jung 2005)
#' @return  A list with the following components:
#' \item{n}{sample size estimate}
#' \item{computed.avepow}{average power}
#' \item{desired.avepow}{desired average power}
#' \item{desired.fdr}{desired FDR}
#' \item{input.pi0}{proportion of tests with a true null hypothesis}
#' \item{alpha}{fixed p-value threshold for multiple testing procedure}
#' \item{n.its}{number of iteration}
#' \item{max.its}{maximum number of iteration, default is 50}
#' \item{n0}{lower limit for initial sample size range}
#' \item{n1}{upper limit for initial sample size range}
#' @references Noether, Gottfried E (1987) Sample size determination for some common nonparametric tests. Journal of the American Statistical Association, 82:645-647.
#' @examples
#' p = rep(c(0.8,0.5),c(100,900));
#' n.fdr.ranksum(fdr = 0.1, pwr = 0.8, p = p, pi0.hat = "BH")
#' @export
n.fdr.ranksum <- function(fdr, pwr, p, pi0.hat = "BH") {
  pi0 <- mean(p == 0.5)
  a <- alpha.power.fdr(fdr, pwr, pi0, pi0.hat)
  res <- find.sample.size(alpha = a, pwr = pwr, avepow.func = average.power.ranksum, p = p)
  res$desired.fdr=fdr
  res$input.pi0=pi0
  res.names=c("n",
              "computed.avepow",
              "desired.avepow",
              "desired.fdr",
              "input.pi0",
              "alpha",
              "n0","n1",
              "n.its",
              "max.its")
  res=res[res.names]
  return(res)
}


#' @title Sample size calculation for one-way ANOVA
#' @description
#' Find the sample size needed to have a desired false discovery rate and average power for a large number of one-way ANOVA tests.
#' @param fdr desired FDR (scalar numeric)
#' @param pwr desired average power (scalar numeric)
#' @param theta sum of ((group mean - overall mean)/stdev)^2 across all groups for each hypothesis test (vector)
#' @param pi0.hat method to estimate proportion \code{pi0} of tests with true null, including: "HH" (p-value histogram height), "HM" (p-value histogram mean), "BH" (Benjamini & Hochberg 1995), "Jung" (Jung 2005)
#' @param k the number of groups to be compared
#' @return  A list with the following components:
#' \item{n}{per-group sample size estimate}
#' \item{computed.avepow}{average power}
#' \item{desired.avepow}{desired average power}
#' \item{desired.fdr}{desired FDR}
#' \item{input.pi0}{proportion of tests with a true null hypothesis}
#' \item{alpha}{fixed p-value threshold for multiple testing procedure}
#' \item{n.its}{number of iteration}
#' \item{max.its}{maximum number of iteration, default is 50}
#' \item{n0}{lower limit for initial sample size range}
#' \item{n1}{upper limit for initial sample size range}
#' @examples
#' theta=rep(c(2,0),c(100,900));
#' n.fdr.oneway(fdr = 0.1, pwr = 0.8, theta = theta, k = 2, pi0.hat = "BH")
#' @export
n.fdr.oneway <- function(fdr, pwr, theta, k, pi0.hat = "BH") {
  pi0 <- mean(theta == 0)
  a <- alpha.power.fdr(fdr, pwr, pi0, pi0.hat)
  res <- find.sample.size(alpha = a, pwr = pwr, avepow.func = average.power.oneway, theta = theta, k = k)
  res$desired.fdr=fdr
  res$input.pi0=pi0
  res.names=c("n",
              "computed.avepow",
              "desired.avepow",
              "desired.fdr",
              "input.pi0",
              "alpha",
              "n0","n1",
              "n.its",
              "max.its")
  res=res[res.names]
  return(res)
}


#' @title Sample size calculation for the Cox proportional hazards regression model
#' @description
#' Find number of events needed to have a desired false discovery rate and average power for a large number of Cox regression models with non-binary covariates.
#' @param fdr desired FDR (scalar numeric)
#' @param pwr desired average power (scalar numeric)
#' @param logHR log hazard ratio (vector)
#' @param v variance of predictor variable (vector)
#' @param pi0.hat method to estimate proportion \code{pi0} of tests with true null, including: "HH" (p-value histogram height), "HM" (p-value histogram mean), "BH" (Benjamini & Hochberg 1995), "Jung" (Jung 2005)
#' @return  A list with the following components:
#' \item{n}{number of events estimate}
#' \item{computed.avepow}{average power}
#' \item{desired.avepow}{desired average power}
#' \item{desired.fdr}{desired FDR}
#' \item{input.pi0}{proportion of tests with a true null hypothesis}
#' \item{alpha}{fixed p-value threshold for multiple testing procedure}
#' \item{n.its}{number of iteration}
#' \item{max.its}{maximum number of iteration, default is 50}
#' \item{n0}{lower limit for initial sample size range}
#' \item{n1}{upper limit for initial sample size range}
#' @note For the test with power calculation based on asymptotic normal approximation, we suggest checking \code{FDRsamplesize2} calculation by simulation.
#' @examples
#' log.HR=log(rep(c(1,2),c(900,100)))
#' v=rep(1,1000)
#' n.fdr.coxph(fdr=0.1, pwr=0.8,logHR=log.HR, v=v, pi0.hat="BH")
#' @references Hsieh, FY and Lavori, Philip W (2000) Sample-size calculations for the Cox proportional hazards regression model with non-binary covariates. Controlled Clinical Trials 21(6):552-560.
#' @export
n.fdr.coxph <- function(fdr, pwr, logHR, v, pi0.hat = "BH") {
  pi0 <- mean(logHR == 0)
  a <- alpha.power.fdr(fdr, pwr, pi0, pi0.hat)
  res <- find.sample.size(a, pwr, average.power.coxph, logHR = logHR, v = v)
  res$desired.fdr=fdr
  res$input.pi0=pi0
  res.names=c("n",
              "computed.avepow",
              "desired.avepow",
              "desired.fdr",
              "input.pi0",
              "alpha",
              "n0","n1",
              "n.its",
              "max.its")
  res=res[res.names]
  return(res)
}


#' @title Sample size calculation for Negative Binomial data
#' @description
#' Find the sample size needed to have a desired false discovery rate and average power for a large number of Negative Binomial comparisons.
#' @param fdr desired FDR (scalar numeric)
#' @param pwr desired average power (scalar numeric)
#' @param log.fc log fold-change (vector), usual null hypothesis is log.fc=0
#' @param mu read depth per gene (vector, same length as log.fc)
#' @param sig coefficient of variation (CV) per gene (vector, same length as log.fc)
#' @param pi0.hat method to estimate proportion \code{pi0} of tests with true null, including: "HH" (p-value histogram height), "HM" (p-value histogram mean), "BH" (Benjamini & Hochberg 1995), "Jung" (Jung 2005)
#' @references SN Hart, TM Therneau, Y Zhang, GA Poland, and J-P Kocher (2013). Calculating Sample Size Estimates for RNA Sequencing Data. Journal of Computational Biology 20: 970-978.
#' @return  A list with the following components:
#' \item{n}{per-group sample size estimate}
#' \item{computed.avepow}{average power}
#' \item{desired.avepow}{desired average power}
#' \item{desired.fdr}{desired FDR}
#' \item{input.pi0}{proportion of tests with a true null hypothesis}
#' \item{alpha}{fixed p-value threshold for multiple testing procedure}
#' \item{n.its}{number of iteration}
#' \item{max.its}{maximum number of iteration, default is 50}
#' \item{n0}{lower limit for initial sample size range}
#' \item{n1}{upper limit for initial sample size range}
#' @note For the test with power calculation based on asymptotic normal approximation, we suggest checking \code{FDRsamplesize2} calculation by simulation.
#' @examples
#' logFC = log(rep(c(1,2),c(900,100)));
#' mu = rep(5,1000);
#' sig = rep(0.6,1000);
#' n.fdr.negbin(fdr = 0.1, pwr = 0.8, log.fc = logFC, mu = mu, sig = sig, pi0.hat = "BH")
#' @export
n.fdr.negbin <- function(fdr, pwr, log.fc, mu, sig, pi0.hat = "BH") {
  pi0 <- mean(log.fc == 0)
  a <- alpha.power.fdr(fdr, pwr, pi0, pi0.hat)
  res <- find.sample.size(a, pwr, average.power.hart,
                          log.fc = log.fc, mu = mu, sig = sig
  )
  res$desired.fdr=fdr
  res$input.pi0=pi0
  res.names=c("n",
              "computed.avepow",
              "desired.avepow",
              "desired.fdr",
              "input.pi0",
              "alpha",
              "n0","n1",
              "n.its",
              "max.its")
  res=res[res.names]
  return(res)
}


#' @title Sample size calculation for t-tests
#' @description
#' Find the sample size needed to have a desired false discovery rate and average power for a large number of t-tests.
#' @param fdr desired FDR (scalar numeric)
#' @param pwr desired average power (scalar numeric)
#' @param delta difference of population means (vector)
#' @param sigma standard deviation (vector or scalar)
#' @param type type of t-test
#' @param alternative one- or two-sided test
#' @param pi0.hat method to estimate proportion \code{pi0} of tests with true null, including: "HH" (p-value histogram height), "HM" (p-value histogram mean), "BH" (Benjamini & Hochberg 1995), "Jung" (Jung 2005)
#' @return  A list with the following components:
#' \item{n}{sample size (per group) estimate}
#' \item{computed.avepow}{average power}
#' \item{desired.avepow}{desired average power}
#' \item{desired.fdr}{desired FDR}
#' \item{input.pi0}{proportion of tests with a true null hypothesis}
#' \item{alpha}{fixed p-value threshold for multiple testing procedure}
#' \item{n.its}{number of iteration}
#' \item{max.its}{maximum number of iteration, default is 50}
#' \item{n0}{lower limit for initial sample size range}
#' \item{n1}{upper limit for initial sample size range}
#' @examples
#' d = rep(c(2,0),c(100,900));
#' n.fdr.ttest(fdr = 0.1, pwr = 0.8, delta = d)
#' @export
n.fdr.ttest <- function(fdr, pwr, delta, sigma = 1, type = "two.sample", pi0.hat = "BH", alternative = "two.sided") {
  pi0 <- mean(delta == 0)
  a <- alpha.power.fdr(fdr, pwr, pi0, pi0.hat)
  res <- find.sample.size(a, pwr, average.power.t.test,
                          delta = delta, sigma = sigma, type = type, alternative = alternative
  )
  res$desired.fdr=fdr
  res$input.pi0=pi0
  res.names=c("n",
              "computed.avepow",
              "desired.avepow",
              "desired.fdr",
              "input.pi0",
              "alpha",
              "n0","n1",
              "n.its",
              "max.its")
  res=res[res.names]
  return(res)
}


#' @title Determines the sample size needed to achieve the desired FDR and average power
#' @description
#' Determines the sample size needed to achieve the desired FDR and average power by given the proportion of true null hypothesis.
#' @param alpha the fixed p-value threshold  (scalar numeric)
#' @param pwr desired average power (scalar numeric)
#' @param avepow.func an R function to compute average power
#' @param n0 lower limit for initial sample size range
#' @param n1 upper limit for initial sample size range
#' @param max.its maximum number of iterations
#' @param ... additional arguments to average power function
#' @return  A list with the following components:
#' \item{n}{a sample size estimate}
#' \item{computed.avepow}{average power}
#' \item{desired.avepow}{desired average power}
#' \item{alpha}{fixed p-value threshold for multiple testing procedure}
#' \item{n.its}{number of iteration}
#' \item{max.its}{maximum number of iteration, default is 50}
#' \item{n0}{lower limit for initial sample size range}
#' \item{n1}{upper limit for initial sample size range}
#' @note For the test with power calculation based on asymptotic normal approximation, we suggest checking \code{FDRsamplesize2} calculation by simulation.
#' @examples
#' #Here, calculating the sample size for the study involving many sign tests
#' average.power.signtest;
#' p.adj = 0.001;
#' p = rep(c(0.8,0.5), c(100,9900));
#' find.sample.size(alpha = p.adj, pwr = 0.8, avepow.func = average.power.signtest, p = p)
#' @export
find.sample.size <- function(alpha, pwr, avepow.func, n0 = 3, n1 = 6, max.its = 50, ...) {
  # compute average power for sample size limits n0 and n1
  pwr0 <- avepow.func(n0, alpha, ...)
  pwr1 <- avepow.func(n1, alpha, ...)

  # double sample size limit until desired average power is reached
  n.its <- 0
  while ((pwr1 < pwr) && (n.its < max.its)) {
    pwr0 <- pwr1
    n0 <- n1
    n1 <- 2 * n1
    pwr1 <- avepow.func(n1, alpha, ...)
    n.its <- n.its + 1
  }

  # use bisection on sample size limits to get sample size
  while (((n1 - n0) > 1) && (n.its < max.its)) {
    n.mid <- floor((n1 + n0) / 2)
    pwr.mid <- avepow.func(n.mid, alpha, ...)
    n.its <- n.its + 1
    if (pwr.mid > pwr) {
      n1 <- n.mid
      pwr1 <- pwr.mid
    } else {
      n0 <- n.mid
      pwr0 <- pwr.mid
    }
  }

  res <- list(
    n = n1,
    computed.avepow = pwr1,
    desired.avepow = pwr,
    alpha = alpha,
    n.its = n.its,
    max.its = max.its,
    n0 = n0, n1 = n1
  )
  # other.input=list(...))

  return(res)
}


#' @title Compute p-value threshold for given the proportion pi0 of tests with a true null
#' @description
#' Given the proportion pi0 of tests with a true null, find the p-value threshold that results in a desired FDR and average power.
#' @param fdr desired FDR (scalar numeric)
#' @param pwr desired average power (scalar numeric)
#' @param pi0 the proportion of tests with a true null hypothesis
#' @param method method to estimate proportion \code{pi0} of tests with true null, including: "HH" (p-value histogram height) , "HM" (p-value histogram mean), "BH" (Benjamini & Hochberg 1995), "Jung" (Jung 2005)
#' @return The fixed p-value threshold for multiple testing procedure
#' @details
#' To get the fixed p-value threshold for multiple testing procedure, 4 approximation methods are provided, they are Benjamini & Hochberg procedure (1995), Jung's formula (2005),
#' method of using p-value histogram height (HH) and method of using p-value histogram mean (HM). For last two methods' details, see Ni Y, Onar-Thomas A, Pounds S.  "Computing Power and Sample Size for the False Discovery Rate in Multiple Applications"
#' @references
#' Pounds S and Cheng C, "Sample size determination for the false discovery rate." Bioinformatics 21.23 (2005): 4263-4271.
#'
#' Gadbury GL, et al. (2004) Power and sample size estimation in high dimensional biology. Statistical Methods in Medical Research 13(4):325-38.
#'
#' Jung,Sin-Ho."Sample size for FDR-control in microarray data analysis." Bioinformatics 21.14 (2005): 3097-3104.
#'
#' Ni Y, Seffernick A, Onar-Thomas A, Pounds S.  "Computing Power and Sample Size for the False Discovery Rate in Multiple Applications", Manuscript.
#'
#' @examples
#' alpha.power.fdr(fdr = 0.1, pwr = 0.9, pi0=0.9, method = "HH")
#' @export
alpha.power.fdr <- function(fdr, pwr, pi0, method = "HH") {
  res <- NA

  # 3-part histogram height
  if (method == "HH") {
    a <- -(1 - fdr) * pi0
    b <- (fdr * (1 - pi0) * pwr + (1 - fdr) * pi0 + (1 - pi0) * (1 - pwr))
    c0 <- -(1 - pi0) * pwr * fdr

    res <- (-b + sqrt(b^2 - 4 * a * c0)) / (2 * a)
    return(res)
  }

  # 3-part histogram expected value
  if (method == "HM") {
    a <- (1 - pi0)
    b <- (1 - pwr) * (1 - pi0) + pi0 * (1 - fdr)
    c0 <- -fdr * (1 - pi0) * pwr
    res <- (-b + sqrt(b^2 - 4 * a * c0)) / (2 * a)
    return(res)
  }

  # Benjamini & Hochberg 1995
  if (method == "BH") {
    res <- (1 - pi0) * pwr * fdr / (1 - fdr * pi0)
  }

  # Jung 2005
  if (method == "Jung") {
    res <- (1 - pi0) * pwr * fdr / (pi0 * (1 - fdr))
  }

  if (is.na(res)) {
    stop('No calculation performed; method must be "HH", "HM", "BH", or "Jung".')
  }

  return(res)
}



#' @title Compute FDR and average power for a given sample size and effect size vector
#' @description
#' For a given fixed sample size and effect size vector,compute FDR and average power as a function of the p-value threshold alpha.
#' @param n sample size
#' @param avepow.func function to compute average power
#' @param null.hypo string to evaluate null hypothesis
#' @param alpha p-value threshold(s) to consider
#' @param method method to estimate proportion pi0 of tests with a true null hypothesis, including: "HH" (p-value histogram height) , "HM" (p-value histogram mean), "BH" (Benjamini & Hochberg 1995), "Jung" (Jung 2005)
#' @param ... additional arguments, including effect size vector for average power function
#' @return  A list with the following components:
#' \item{n}{input sample size}
#' \item{avepow.func}{average power function}
#' \item{null.hypo}{null hypothesis string}
#' \item{pi0}{computed value of pi0}
#' \item{method}{method to estimate proportion \code{pi0} of tests with true null, including: "HH" (p-value histogram height) , "HM" (p-value histogram mean), "BH" (Benjamini & Hochberg 1995), "Jung" (Jung 2005)}
#' \item{other.args}{additional arguments}
#' \item{res.tbl}{table of alpha, fdr, and average power}
#' @references
#' Pounds S and Cheng C, "Sample size determination for the false discovery rate." Bioinformatics 21.23 (2005): 4263-4271.
#'
#' Gadbury GL, et al. (2004) Power and sample size estimation in high dimensional biology. Statistical Methods in Medical Research 13(4):325-38.
#'
#' Jung,Sin-Ho."Sample size for FDR-control in microarray data analysis." Bioinformatics 21.14 (2005): 3097-3104.
#'
#' Ni Y, Seffernick A, Onar-Thomas A, Pounds S.  "Computing Power and Sample Size for the False Discovery Rate in Multiple Applications", Manuscript.
#'
#' @examples
#' n = 50; # number of events
#' logHR = rep(c(0,0.5),c(950,50));
#' v = rep(1,length(logHR));  # variance of predictor variable (vector)
#' res = fdr.avepow(n,average.power.coxph,"logHR==0",logHR=logHR,v=v);
#' res$pi0;
#' head(res$res.tbl)
#' @export
fdr.avepow=function(n,                 # sample size
                    avepow.func,       # function to compute average power
                    null.hypo,         # string to evaluate null hypothesis
                    alpha=1:100/1000,  # p-value threshold(s) to consider
                    method="BH",       # method to estimate proportion $\pi_0$ of tests with true null
                    ...)               # additional named arguments, including effect size vector for average power function

{
  alpha=alpha[(alpha>0)&(alpha<1)]
  if (length(alpha)<1)
    stop("Invalid alpha.")

  n.alpha=length(alpha)
  eff.size=eval(parse(text=null.hypo))
  pi0=mean(eff.size)
  if (is.na(pi0))
    stop("Invalid null.hypo string.")

  pwr=rep(NA,length(alpha))
  for (i in 1:n.alpha)  pwr[i]=avepow.func(n,alpha[i],...)
  fdr=fdr.power.alpha(alpha,pwr,pi0,method)

  res.tbl=cbind(alpha=alpha,fdr=fdr,pwr=pwr)

  res=list(n=n,                      # input sample size
           avepow.func=avepow.func,  # average power function
           null.hypo=null.hypo,      # null hypothesis string
           pi0=pi0,                  # computed value of pi0
           method=method,            # method to estimate $pi_0$
           other.args=list(...),     # additional arguments
           res.tbl=res.tbl)          # table of alpha, fdr, and average power

  return(res)
}


#' @title Compute FDR for given p-value threshold, average power and proportion of tests with a true null
#' @description
#' Compute the FDR for given values of the p-value threshold alpha, average power, and proportion pi0 of tests with a true null hypothesis.
#' @param alpha p-value threshold (vector)
#' @param pwr average power
#' @param pi0 actual proportion of tests with a true null hypothesis
#' @param method method to estimate proportion \code{pi0} of tests with true null, including: "HH" (p-value histogram height), "HM" (p-value histogram mean), "BH" (Benjamini & Hochberg 1995), "Jung" (Jung 2005)
#' @return  FDR
#' @references
#' Pounds S and Cheng C, "Sample size determination for the false discovery rate." Bioinformatics 21.23 (2005): 4263-4271.
#'
#' Gadbury GL, et al. (2004) Power and sample size estimation in high dimensional biology. Statistical Methods in Medical Research 13(4):325-38.
#'
#' Jung,Sin-Ho."Sample size for FDR-control in microarray data analysis." Bioinformatics 21.14 (2005): 3097-3104.
#'
#' Ni Y, Seffernick A, Onar-Thomas A, Pounds S.  "Computing Power and Sample Size for the False Discovery Rate in Multiple Applications", Manuscript.
#'
#' @examples
#' alpha = 1:100/1000;
#' pwr = rep(0.8,length(alpha));
#' pi0 = 0.95;
#' fdr.power.alpha(alpha,pwr,pi0,method="HH")
#' @importFrom stats weighted.mean
#' @export
fdr.power.alpha=function(alpha,        # p-value threshold (vector)
                         pwr,          # average power
                         pi0,          # actual proportion of tests with a true null hypothesis
                         method="HH")  # method to estimate the proportion of tests with a true null hypothesis
{
  if (method=="BH")
  {
    fdr=alpha/(pi0*alpha+(1-pi0)*pwr)
    return(fdr)
  }

  if (method=="HH")
  {
    pi0.hat=pi0+(1-pi0)*(1-pwr)/(1-alpha)
    fdr=pi0.hat*alpha/(pi0*alpha+(1-pi0)*pwr)
    return(fdr)
  }

  if (method=="HM")
  {
    mns=c(alpha/2,(1+alpha)/2)
    wgt=c(pi0*alpha+(1-pi0)*pwr,
          pi0*(1-alpha)+(1-pi0)*(1-pwr))
    mn.pval=weighted.mean(mns,wgt)
    pi0.hat=2*mn.pval
    fdr=pi0.hat*alpha/(pi0*alpha+(1-pi0)*pwr)
    return(fdr)
  }

  if (method=="Jung")
  {
    pi0.hat=pi0
    fdr=pi0.hat*alpha/(pi0*alpha+(1-pi0)*pwr)
    return(fdr)
  }

  stop('Invalid pi0 method request: method must be "BH", "HM", "HH", or "Jung".')

}

