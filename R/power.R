#' @title Compute the power of a single-predictor Cox regression model
#' @description Use the formula of Hsieh and Lavori (2000) to compute the power of a single-predictor Cox model
#' @param n number of events (scalar)
#' @param alpha p-value threshold (scalar)
#' @param logHR log hazard ratio (vector)
#' @param v variance of predictor variable (vector)
#' @return Vector of power estimates for two-sided test
#' @references Hsieh, FY and Lavori, Philip W (2000) Sample-size calculations for the Cox proportional hazards regression model with non-binary covariates. Controlled Clinical Trials 21(6):552-560.
#' @examples
#' logHR = log(rep(c(1, 2),c(900, 100)));
#' v = rep(1, 1000);
#' res = power.cox(n = 50,alpha = 0.05,logHR = logHR, v = v)
#' @importFrom stats pnorm qnorm
#' @export
power.cox <- function(n,
                      alpha,
                      logHR,
                      v) {
  pnorm(qnorm(alpha / 2), sqrt(n * v) * logHR) +
    1 - pnorm(qnorm(1 - alpha / 2), sqrt(n * v) * logHR)
}


#' @title Compute power for RNA-seq experiments assuming Negative Binomial distribution
#' @description
#' Use the formula of Hart et al (2013) to compute power for comparing RNA-seq expression across two groups assuming a Negative Binomial distribution
#' @param n per-group sample size (scalar)
#' @param alpha p-value threshold (scalar)
#' @param log.fc log fold-change (vector), usual null hypothesis is log.fc=0
#' @param mu read depth per gene (vector, same length as log.fc)
#' @param sig coefficient of variation (CV) per gene (vector, same length as log.fc)
#' @details
#' This function is based on equation (1) of Hart et al (2013). It assumes a Negative Binomial model for RNA-seq read counts and equal sample size per group.
#' @return Vector of power estimates for the set of two-sided tests
#' @references SN Hart, TM Therneau, Y Zhang, GA Poland, and J-P Kocher (2013). Calculating Sample Size Estimates for RNA Sequencing Data. Journal of Computational Biology 20: 970-978.
#' @examples
#' n.hart = 2*(qnorm(0.975)+qnorm(0.9))^2*(1/20+0.6^2)/(log(2)^2)   # Equation (6) of Hart et al
#' power.hart(n.hart,0.05,log(2),20,0.6)                            # Recapitulate 90% power
#' @importFrom stats pnorm qnorm
#' @export
power.hart <- function(n,
                       alpha,
                       log.fc,
                       mu,
                       sig) {
  z.alpha <- qnorm(alpha / 2)
  res <- pnorm(z.alpha, sqrt(n * log.fc^2 / (2 * (1 / mu + sig^2)))) +
    pnorm(-z.alpha, sqrt(n * log.fc^2 / (2 * (1 / mu + sig^2))), lower.tail = F)
  return(res)
}


#' @title Compute power of the rank-sum test
#' @description
#' Compute power of rank-sum test; Uses formula of Noether (JASA 1987)
#' @param n per-group sample size (scalar)
#' @param alpha p-value threshold (scalar)
#' @param p Pr(Y>X), as in Noether (JASA 1987)
#' @details
#' In most applications, the null effect size will be designated by p = 0.5
#' @return Vector of power estimates for two-sided tests
#' @references Noether, Gottfried E (1987) Sample size determination for some common nonparametric tests. Journal of the American Statistical Association, 82:645-647.
#' @examples
#' p = rep(c(0.8,0.5),c(100,900))
#' res = power.ranksum(n = 50, alpha = 0.5, p=p)
#' @importFrom stats qnorm pnorm
#' @export
power.ranksum <- function(n, alpha, p) {
  mu0 <- 0.5 * n * n
  mu1 <- p * n * n
  delta <- (mu1 - mu0)
  sig2 <- n * n * (2 * n + 1) / 12
  sig <- sqrt(sig2)
  z.reject <- -abs(qnorm(alpha / 2, 0, sig))
  pow <- pnorm(z.reject, delta, sig) + pnorm(-z.reject, delta, sig, lower.tail = F)
  return(pow)
}


#' @title Compute power of one-way ANOVA
#' @description
#' Compute power of one-way ANOVA; Uses classical power formula for ANOVA; Assumes equal variance and sample size
#' @param n per-group sample size (scalar)
#' @param alpha p-value threshold (scalar)
#' @param theta sum of ((group mean - overall mean)/stdev)^2 across all groups for each hypothesis test(vector)
#' @param k the number of groups to be compared, default k=2
#' @details
#' For many applications, the null effect is zero for the parameter theta described above
#' @returns Vector of power estimates for test of equal means
#' @examples
#' theta=rep(c(2,0),c(100,900));
#' res = power.oneway(n = 50, alpha = 0.05, theta = theta, k = 2)
#' @importFrom stats pf qf
#' @export
power.oneway <- function(n, alpha, theta, k = 2) {
  1 - pf(qf(1 - alpha, k - 1, k * (n - 1)), k - 1, k * (n - 1), ncp = n * theta)
}


#' @title Compute power of the sign test
#' @description
#' Use the Noether (1987) formula to compute the power of the sign test
#' @param n sample size (scalar)
#' @param alpha p-value threshold (scalar)
#' @param p Pr(X>0), as in Noether (JASA 1987)
#' @details
#' In most applications, the null effect size will be designated by p = 0.5
#' @return Vector of power estimates for two-sided tests
#' @references Noether, Gottfried E (1987) Sample size determination for some common nonparametric tests. Journal of the American Statistical Association, 82:645-647.
#' @examples
#' p = rep(c(0.8,0.5),c(100,900));
#' res = power.signtest(n = 50, alpha = 0.05, p = p)
#' @export
power.signtest <- function(n, alpha, p) {
  mu0 <- 0.5 * n
  mu1 <- p * n
  sig0 <- sqrt(n * 0.25)
  sig1 <- sqrt(n * p * (1 - p))
  pow <- pnorm(qnorm(alpha / 2, mu0, sig0), mu1, sig1) +
    pnorm(qnorm(1 - alpha / 2, mu0, sig0), mu1, sig1, lower.tail = F)
  return(pow)
}


#' @title Compute power of the signed-rank test
#' @description
#' Use the Noether (1987) formula to compute the power of the signed-rank test
#' @param n sample size (scalar)
#' @param alpha p-value threshold (scalar)
#' @param p1 Pr(X>0), as in Noether (JASA 1987)
#' @param p2 Pr(X+X'>0), as in Noether (JASA 1987)
#' @details
#' In most applications, the null effect size will be designated by p1 = p2 = 0.5
#' @return Vector of power estimates for two-sided tests
#' @references Noether, Gottfried E (1987) Sample size determination for some common nonparametric tests. Journal of the American Statistical Association, 82:645-647.
#' @examples
#' p1 = rep(c(0.8,0.5),c(100,900));
#' p2 = rep(c(0.8,0.5),c(100,900));
#' res = power.signrank(n = 50, alpha = 0.05, p1 = p1, p2 = p2)
#' @importFrom stats pnorm qnorm
#' @export
power.signrank <- function(n,
                           alpha,
                           p1,
                           p2) {
  mu0 <- 0.25 * n * (n + 1)
  mu1 <- n*p1 + 0.5 * p2 * n * (n - 1)
  delta <- (mu1 - mu0)
  sig2 <- n * (n + 1) * (2 * n + 1) / 24
  sig <- sqrt(sig2)
  z.reject <- -abs(qnorm(alpha / 2, 0, sig))
  pow <- pnorm(z.reject, delta, sig) + pnorm(-z.reject, delta, sig, lower.tail = F)
  return(pow)
}


#' @title Compute power for RNA-Seq experiments assuming Poisson distribution
#' @description
#' Use the formula of Li et al (2013) to compute power for comparing RNA-seq expression across two groups assuming the Poisson distribution
#' @param n per-group sample size
#' @param alpha p-value threshold (scalar)
#' @param rho fold-change, usual null hypothesis is that rho=1 (vector)
#' @param mu0 average count in control group
#' @param w ratio of the total number of reads mapped between the two groups
#' @param type type of test: "w" for Wald, "s" for score, "lw" for log-transformed Wald, "ls" for log-transformed score
#' @details
#' This function computes the power for each of a series of two-sided tests defined by the input parameters. The power is based on the sample size formulas in equations (10-13) of Li et al (2013). Also, note that the null.effect is set to 1 in the examples because the usual null hypothesis is that the fold-change = 1.
#' @return Vector of power estimates for two-sided tests
#' @references C-I Li, P-F Su, Y Guo, and Y Shyr (2013). Sample size calculation for differential expression analysis of RNA-seq data under Poisson distribution. Int J Comput Biol Drug Des 6(4). <doi:10.1504/IJCBDD.2013.056830>
#' @examples
#' power.li(n = 88, alpha = 0.05, rho = 1.25, mu0 = 5, w = 0.5,type = "w")
#' # recapitulate 80% power in Table 1 of Li et al (2013)
#' @importFrom stats pnorm qnorm
#' @export
power.li <- function(n,
                     alpha,
                     rho,
                     mu0,
                     w ,
                     type)
{
  z.alpha <- qnorm(alpha / 2)
  if (!is.element(type, c("w", "lw", "s", "ls"))) {
    stop("type must be 'w', 's', 'ls', or 'lw'.")
  }
  if (any(w <= 0)) stop("w must be >0.")
  if (any(rho <= 0)) stop("rho must be >0.")
  if (type == "w") {
    mu.z <- sqrt(n * mu0 * (rho - 1)^2 / (1 + rho / w))
    res <- pnorm(z.alpha, mu.z) + pnorm(-z.alpha, mu.z, lower.tail = F)
    return(res)
  }
  if (type == "lw") {
    mu.z <- sqrt(n * mu0 * (log(rho)^2) / (1 + 1 / (rho * w)))
    res <- pnorm(z.alpha, mu.z) + pnorm(-z.alpha, mu.z, lower.tail = F)
    return(res)
  }
  if (type == "s") {
    z.cut <- z.alpha * sqrt((1 + w * rho) / (w + rho))
    mu.z <- sqrt(n * mu0 * (rho - 1)^2 / (1 + rho / w))
    res <- pnorm(z.cut, mu.z) + pnorm(-z.cut, mu.z, lower.tail = F)
    return(res)
  }
  if (type == "ls") {
    z.cut <- z.alpha * sqrt((2 + w + 1 / w) / ((1 + w * rho) * (1 + 1 / (w * rho))))
    mu.z <- sqrt(n * mu0 * (log(rho)^2) / (1 + 1 / (rho * w)))
    res <- pnorm(z.cut, mu.z) + pnorm(-z.cut, mu.z, lower.tail = F)
    return(res)
  }
}


#' @title Compute power for Fisher's exact test
#' @param p1  probability in one group (scalar)
#' @param p2  probability in other group (scalar)
#' @param n per-group sample size (scalar)
#' @param alpha p-value threshold (scalar)
#' @param alternative one- or two-sided test, must be one of "greater", "less", or "two.sided"
#' @return Power estimate for one- or two-sided tests
#' @examples
#' power.fisher(p1 = 0.5, p2 = 0.9, n=20, alpha = 0.05, alternative = 'two.sided')
#' @export
#' @importFrom stats dbinom fisher.test
power.fisher <- function(p1, p2, n, alpha, alternative)
{
  pr1 <- dbinom(0:n, n, p1, log = T)
  pr2 <- dbinom(0:n, n, p2, log = T)
  y1 <- rep(0:n, each = n + 1)
  y2 <- rep(0:n, times = n + 1)
  y <- cbind(y1, n - y1, y2, n - y2)
  pr <- exp(rep(pr1, each = n + 1) + rep(pr2, times = n + 1))
  res <- 0
  for (i in 1:length(pr))
  {
    pval <- fisher.test(matrix(y[i, ], 2, 2), alternative = alternative)$p.value
    res <- res + (pval <= alpha) * pr[i]
  }

  return(res)
}


#' @title Compute power of the t-test for non-zero correlation
#' @param n sample size (scalar)
#' @param alpha p-value threshold (scalar)
#' @param rho population correlation coefficient (vector)
#' @return Vector of power estimates for two-sided tests
#' @details
#' For many applications, the null.effect is rho = 0
#' @examples
#' rho = rep(c(0.3,0),c(100,900));
#' res = power.tcorr(n = 50, alpha = 0.05, rho = rho)
#' @importFrom stats pf qf
#' @export
power.tcorr=function(n,
                     alpha,
                     rho)
{
  1-pf(qf(1-alpha,1,n-2),1,n-2,ncp=(n-2)*rho^2/(1-rho^2))
}
