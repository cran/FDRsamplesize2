
<!-- README.md is generated from README.Rmd. Please edit that file -->

# FDRsamplesize2

<!-- badges: start -->
<!-- badges: end -->

The goal of *FDRsamplesize2* is to computer the average power and
determine the sample size needed to achieve desired FDR and desired
average power. For power calculation, the package covers various
statistical test settings,including:

- two-samples $t$-test
- one-sample $t$-test
- rank-sum test
- signed-rank test
- Fisher’s exact test
- the $t$-test for correlation
- comparison of two Poisson distributions
- comparison of two negative binomial distributions
- the two-proportions $z$-test
- one-way ANOVA
- Cox proportional hazards regression

## Installation

You can install *FDRsamplesize2* from github:

``` r
# install.packages("devtools")
devtools::install_github("yonghui-ni/FDRsamplesize_2.0")
```

Or from CRAN:

``` r
# install.packages("FDEsamplesize2")
```

## Load package

load the package

``` r
library("FDRsamplesize2")
```

## Example

This is an example which shows you how to determine the sample size
necessary to identify genes with microarray expression signals that are
associated with progression-free survival. A single-predictor Cox
regression model will be used to test the association of each gene’s
expression with progression-free survival (PFS). We are interested in
determining the number of events necessary to identify 80% of genes
truly associate with PFS while controlling the FDR at 10% in a setting
in which 1% of the genes have a hazard ratio of 2 and a variance of 1.

``` r
log.HR=log(rep(c(1,2),c(9900,100)))   # log hazard ratio for each gene
v=rep(1,10000)                        # variance of each gene
res=n.fdr.coxph(fdr=0.1, pwr=0.8,logHR=log.HR, v=v, pi0.hat="BH")
res
#> $n
#> [1] 37
#> 
#> $computed.avepow
#> [1] 0.8139159
#> 
#> $desired.avepow
#> [1] 0.8
#> 
#> $desired.fdr
#> [1] 0.1
#> 
#> $input.pi0
#> [1] 0.99
#> 
#> $alpha
#> [1] 0.0008879023
#> 
#> $n0
#> [1] 36
#> 
#> $n1
#> [1] 37
#> 
#> $n.its
#> [1] 7
#> 
#> $max.its
#> [1] 50
```

Step by step to calculate power and sample size. This procedure is more
flexible for user to plug in average power function of other statistical
tests that are not available in *FDRsamplesize2*.

``` r
pi0=9900/10000                      # proportion of true null hypothesis
adj.p=alpha.power.fdr(fdr=0.1,
                       pwr=0.8,
                       pi0=pi0,
                       method="BH")
print(adj.p)
#> [1] 0.0008879023

find.sample.size(alpha=adj.p,                       # Fixed p-value threshold     
                 pwr=0.8,                           # Power of single α-level test
                 avepow.func=average.power.coxph,   # an R function to compute average power  
                 n0=3,                              # Lower limit for initial sample size range
                 n1=6,                              # Upper limit for initial sample size range
                 logHR=log(rep(c(1,2),c(9900,100))),
                 v=rep(1,10000), 
                 max.its=50                         # Maximum number of iterations, default is 50
                 )
#> $n
#> [1] 37
#> 
#> $computed.avepow
#> [1] 0.8139159
#> 
#> $desired.avepow
#> [1] 0.8
#> 
#> $alpha
#> [1] 0.0008879023
#> 
#> $n.its
#> [1] 7
#> 
#> $max.its
#> [1] 50
#> 
#> $n0
#> [1] 36
#> 
#> $n1
#> [1] 37
```

To get the fixed p-value threshold $\alpha$ in *alpha.power.fdr*
function for multiple testing procedure, approximation method options
include:

- Benjamini and Hochberg 1995 (BH)
- Jung 2005 (Jung)
- p-value histogram height (HH)
- p-value histogram mean (HM)
