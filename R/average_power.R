#' @title Compute the average power of many Cox regression models
#' @description
#' Compute the average power of many Cox regression models for a given number of events, p-value threshold, vector of effect sizes (log hazard ratio),and variance of predictor variables
#' @param n number of events (scalar)
#' @param alpha p-value threshold (scalar)
#' @param logHR log hazard ratio (vector)
#' @param v variance of predictor variable (vector)
#' @return Average power estimate for multiple testing procedure
#' @examples
#' logHR = log(rep(c(1, 2),c(900, 100)));
#' v = rep(1, 1000);
#' average.power.coxph(n = 50, alpha = 0.05, logHR = logHR, v = v)
#' @references Hsieh, FY and Lavori, Philip W (2000) Sample-size calculations for the Cox proportional hazards regression model with non-binary covariates. Controlled Clinical Trials 21(6):552-560.
#' @seealso \code{\link{power.cox}} for more details about power calculation of single-predictor Cox regression model. The power calculation is based on asymptotic normal approximation.
#' @export
average.power.coxph <- function(n, alpha, logHR, v) {
  m <- length(logHR)
  if (length(v) == 1) v <- rep(v, m)
  alt <- (logHR != 0)
  logHR <- logHR[alt]
  v <- v[alt]
  res <- mean(power.cox(n, alpha, logHR, v))
  return(res)
}


#' @title Compute average power for RNA-seq experiments assuming Negative Binomial distribution
#' @param n per-group sample size (scalar)
#' @param alpha p-value threshold (scalar)
#' @param log.fc log fold-change (vector), usual null hypothesis is log.fc=0
#' @param mu read depth per gene (vector, same length as log.fc)
#' @param sig coefficient of variation (CV) per gene (vector, same length as log.fc)
#' @return Average power estimate for multiple testing procedure
#' @details
#' The power function is based on equation (1) of Hart et al (2013). It assumes a Negative Binomial model for RNA-seq read counts and equal sample size per group.
#' @references SN Hart, TM Therneau, Y Zhang, GA Poland, and J-P Kocher (2013). Calculating Sample Size Estimates for RNA Sequencing Data. Journal of Computational Biology 20: 970-978.
#' @examples
#' logFC = log(rep(c(1,2),c(900,100)));
#' mu = rep(5,1000);
#' sig = rep(0.6,1000);
#' average.power.hart(n = 50, alpha = 0.05,log.fc = logFC, mu = mu, sig = sig)
#' @seealso \code{\link{power.hart}} for more details about power calculation of data under Negative Binomial distribution. The power calculation is based on asymptotic normal approximation.
#' @export
average.power.hart <- function(n, alpha, log.fc, mu, sig) {
  alt <- (log.fc != 0)
  log.fc <- log.fc[alt]
  pwr <- power.hart(n = n, alpha = alpha, log.fc = log.fc, mu = mu, sig = sig)
  res <- mean(pwr)
  return(res)
}


#' @title Compute average power of rank-sum tests
#' @description
#' Compute average power of rank-sum tests
#' @param n sample size (scalar)
#' @param alpha p-value threshold (scalar)
#' @param p Pr(Y>X), as in Noether (JASA 1987)
#' @return Average power estimate for multiple testing procedure
#' @examples
#' p = rep(c(0.8,0.5),c(100,900));
#' average.power.ranksum(n = 50, alpha = 0.05, p=p)
#' @seealso \code{\link{power.ranksum}} for more details about power calculation of rank-sum test. The power calculation is based on asymptotic normal approximation.
#' @export
average.power.ranksum <- function(n, alpha, p) {
  alt <- (p != 0.5)
  p <- p[alt]
  pwr <- power.ranksum(n, alpha, p)
  res <- mean(pwr)
  return(res)
}


#' @title Compute average power of many one-way ANOVA tests
#' @description
#' Compute average power of many one-way ANOVA tests
#' @param n per-group sample size (scalar)
#' @param alpha p-value threshold (scalar)
#' @param theta sum of ((group mean - overall mean)/stdev)^2 across all groups for each hypothesis test(vector)
#' @param k the number of groups to be compared
#' @return Average power estimate for multiple testing procedure
#' @examples
#' theta=rep(c(2,0),c(100,900));
#' average.power.oneway(n = 50, alpha = 0.05, theta = theta, k = 2)
#' @seealso \code{\link{power.oneway}} for more details about power calculation of one-way ANOVA
#' @export
average.power.oneway <- function(n, alpha, theta, k) {
  alt <- (theta != 0)
  theta <- theta[alt]
  pwr <- power.oneway(n, alpha, theta, k)
  res <- mean(pwr)
  return(res)
}


#' @title Compute average power of many sign tests
#' @param n sample size (scalar)
#' @param alpha p-value threshold (scalar)
#' @param p Pr(Y>X), as in Noether (JASA 1987)
#' @return Average power estimate for multiple testing procedure
#' @examples
#' p = rep(c(0.8,0.5),c(100,900));
#' average.power.signtest(n = 50, alpha = 0.05, p=p)
#' @seealso \code{\link{power.signtest}} for more details about power calculation of sign test. The power calculation is based on asymptotic normal approximation.
#' @importFrom stats qnorm pnorm
#' @export
average.power.signtest <- function(n, alpha, p) {
  alt <- (p != 0.5)
  p <- p[alt]
  pwr <- power.signtest(n, alpha, p)
  res <- mean(pwr)
  return(res)
}


#' @title Compute average power of many signed-rank tests
#' @param n sample size (scalar)
#' @param alpha p-value threshold (scalar)
#' @param p1 Pr(X>0), as in Noether (JASA 1987)
#' @param p2 Pr(X+X'>0), as in Noether (JASA 1987)
#' @return Average power estimate for multiple testing procedure
#' @examples
#' p1 = rep(c(0.8,0.5),c(100,900));
#' p2 = rep(c(0.8,0.5),c(100,900));
#' average.power.signrank(n = 50, alpha = 0.05, p1 = p1, p2 = p2)
#' @seealso \code{\link{power.signrank}} for more details about power calculation of signed-rank test. The power calculation is based on asymptotic normal approximation.
#' @export
average.power.signrank <- function(n, alpha, p1, p2) {
  alt <- intersect(which(p1 == 0.5),which(p2 == 0.5))
  p1 <- p1[-alt]
  p2 <- p2[-alt]
  pwr <- power.signrank(n, alpha=alpha, p1=p1 ,p2=p2)
  res <- mean(pwr)
  return(res)
}



#' @title Compute average power for RNA-Seq experiments assuming Poisson distribution
#' @description
#' Use the formula of Li et al (2013) to compute power for comparing RNA-seq expression across two groups assuming the Poisson distribution.
#' @param n per-group sample size
#' @param alpha p-value threshold (scalar)
#' @param rho fold-change, usual null hypothesis is that rho=1 (vector)
#' @param mu0 average count in control group (vector)
#' @param w ratio of the total number of reads mapped between the two groups (scalar or vector)
#' @param type type of test: "w" for Wald, "s" for score, "lw" for log-transformed Wald, "ls" for log-transformed score
#' @details
#' This function computes the average power for a series of two-sided tests defined by the input parameters. The power is based on the sample size formulas in equations (10-13) of Li et al (2013). Also, note that the null.effect is set to 1 in the examples because the usual null hypothesis is that the fold-change = 1.
#' @return Average power estimate for multiple testing procedure
#' @references C-I Li, P-F Su, Y Guo, and Y Shyr (2013). Sample size calculation for differential expression analysis of RNA-seq data under Poisson distribution. Int J Comput Biol Drug Des 6(4).<doi:10.1504/IJCBDD.2013.056830>
#' @examples
#' rho = rep(c(1,1.25),c(900,100));
#' mu0 = rep(5,1000);
#' w = rep(0.5,1000);
#' average.power.li(n = 50, alpha = 0.05, rho = rho, mu0 = mu0, w = w, type = "w")
#' @seealso \code{\link{power.li}} for more details about power calculation of data under Poisson distribution
#' @export
average.power.li <- function(n, alpha, rho, mu0, w, type) {
  alt <- (rho != 1)
  rho <- rho[alt]
  pwr <- power.li(n, alpha = alpha, rho = rho, mu0 = mu0, w = w, type = type)
  res <- mean(pwr)
  return(res)
}


#' @title Compute average power of many Fisher's exact tests
#' @param p1  probability in one group (vector)
#' @param p2  probability in other group (vector)
#' @param n per-group sample size
#' @param alpha p-value threshold
#' @param alternative one- or two-sided test
#' @return Average power estimate for multiple testing procedure
#' @examples
#' set.seed(1234);
#' p1 = sample(seq(0,0.5,0.1),5,replace = TRUE);
#' p2 = sample(seq(0.5,1,0.1),5,replace = TRUE);
#' average.power.fisher(p1 = p1,p2 = p2,n = 20,alpha = 0.05,alternative = "two.sided")
#' @seealso \code{\link{power.fisher}} for more details about power calculation of Fisher's exact test
#' @importFrom stats dbinom fisher.test
#' @export
average.power.fisher <- function(p1, p2, n, alpha, alternative)
{
  alt <- (p1 != p2)
  p1 <- p1[alt]
  p2 <- p2[alt]
  n <- rep(n, length(p1))

  pwr <- 0
  for (i in 1:length(p1))
  {
    temp <- power.fisher(p1 = p1[i], p2 = p2[i], n = n[i], alpha = alpha, alternative = alternative)

    pwr <- pwr + temp
  }
  pwr <- pwr / length(p1)
  return(pwr)
}


#' @title Compute average power of many t-tests for non-zero correlation
#' @param n sample size (scalar)
#' @param alpha p-value threshold (scalar)
#' @param rho population correlation coefficient (vector)
#' @return Average power estimate for multiple testing procedure
#' @details
#' For many applications, the null.effect is rho = 0
#' @examples
#' rho = rep(c(0.3,0),c(100,900));
#' average.power.tcorr(n = 50, alpha = 0.05, rho = rho)
#' @seealso \code{\link{power.tcorr}} for more details about power calculation of t-test for non-zero correlation
#' @export
average.power.tcorr=function(n,alpha,rho)

{
  alt=(rho!=0)
  rho=rho[alt]
  pwr=power.tcorr(n,alpha,rho)
  res=mean(pwr)
  return(res)
}


#' @title Computer average power of many two proportion z-tests
#' @description
#' Computer average power of many two proportion z-tests.The power calculation of two proportion z-test is based on asymptotic normal approximation.
#' @param n  per-group sample size (scalar)
#' @param p1  probability in one group (vector)
#' @param p2  probability in other group (vector)
#' @param alpha p-value threshold (scalar)
#' @param alternative one- or two-sided test
#' @return Average power estimate for multiple testing procedure
#' @examples
#' set.seed(1234);
#' p1 = sample(seq(0,0.5,0.1),40,replace = TRUE);
#' p2 = sample(seq(0.5,1,0.1),40,replace = TRUE);
#' average.power.twoprop(n = 30, alpha = 0.05, p1 = p1,p2 = p2,alternative="two.sided")
#' @importFrom stats power.prop.test
#' @export
average.power.twoprop <- function(n, alpha, p1, p2, alternative) {
  alt <- (p1 != p2)
  p1 <- p1[alt]
  p2 <- p2[alt]
  pwr <- power.prop.test(n = n, p1 = p1, p2 = p2, sig.level = alpha, alternative = alternative)$power
  res <- mean(pwr)
  return(res)
}

#' @title Compute average power of many t-tests
#' @description
#' Compute average power of many t-tests; Uses classical power formula for t-test; Assumes equal variance and sample size
#' @param n per-group sample size (scalar)
#' @param alpha p-value threshold (scalar)
#' @param delta difference of population means (vector)
#' @param sigma standard deviation (vector or scalar, default=1)
#' @param type type of t-test: "two.sample", "one.sample"
#' @param alternative one- or two-sided test
#' @return Average power estimate for multiple testing procedure
#' @examples
#' d = rep(c(2,0),c(100,900));
#' average.power.t.test(n = 20, alpha = 0.05,delta = d)
#' @importFrom stats power.t.test
#' @export
average.power.t.test <- function(n, alpha, delta, sigma=1, type ="two.sample", alternative="two.sided") {
  alt <- (delta != 0)
  delta <- delta[alt]
  pwr <- power.t.test(n=n, delta=delta, sd=sigma, sig.level=alpha, type=type, alternative=alternative)
  res <- mean(pwr$power)
  return(res)
}


