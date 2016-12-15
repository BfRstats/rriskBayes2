# rriskBayes2

## Description

This packages provides a collection of functions for fitting Bayesian models (related to the rrisk project). The functions can be used as stand-alone applications or launched during an rrisk session. 
The following Bayesian models are implemented in this package:

* **PEM** - prevalence estimation under misclassification
* **ZIP** - estimation of a zero inflated Poisson model
* **ZINB** - zero inflated negative binomial model

## Details

Collection of functions for fitting Bayesian models. The functions can be used as stand-alone application or launched during an 'rrisk' session.

This package is a part of the rrisk project and contains functions for fitting Bayesian models using the package JAGS. This package does not depend on the whole rrisk project and can be used separately. The [rrisk project](http://www.bfr.bund.de/cd/52158) has been developed at the BfR (Federal Institute for Risk Assessment) and is available at https://github.com/BfRstats/rrisk.

## Author(s)

* Natalia Belgorodski (at the time STAT-UP Statistical Consulting)
* Matthias Greiner matthias.greiner@bfr.bund.de (Federal Institute for Risk Assessment, Germany)
* Alexander Engelhardt alexander.engelhardt@stat-up.com (STAT-UP Statistical Consulting)
* Christina Yassouridis christina.yassouridis@stat-up.com (STAT-UP Statistical Consulting)

## Installation

Use the following in R: 

```
if (!require(devtools))
{
  install.packages("devtools")
  library(devtools)
}

install_github("BfRstats/rriskBayes2", dependencies = TRUE)
```

## Examples

````
library(rriskBayes2)

#------------------------------------------
# Example of PEM model (k>1)
#------------------------------------------
pi <- 0.01
se <- 0.96
se.n <- 1000
sp <- 0.99
sp.n <- 1000
n <- sample(10:1000,1,replace=TRUE)  # stochatsic sample size
k <- sample(5:50,1,replace=FALSE)    # stochastic pool size

# Parameters for beta priors
se.a <- se.n*se+1
se.b <- se.n*(1-se)+1
sp.a <- sp.n*sp+1
sp.b <- sp.n*(1-sp)+1

# Random number of positive pools (x) considering uncertainty of se and sp
ap <- pi*se + (1-pi)*(1-sp)
p.pos <- 1-(1-ap)^k
x <- rbinom(1,prob=p.pos,size=n)

# Estimate using Bayes model at individual level
resPEM1 <- rrisk.BayesPEM(x=x, n=n,k=k,
     prior.pi=c(1,1),prior.se=c(se.a,se.b),prior.sp=c(sp.a,sp.b),
     misclass="individual")
resPEM1@results

# Estimate using Bayes model at pool level
resPEM2 <- rrisk.BayesPEM(x=x, n=n,k=k,
     prior.pi=c(1,1),prior.se=c(se.a,se.b),prior.sp=c(sp.a,sp.b),
     misclass="pool")
resPEM2@results

# Estimate using Bayes model compared
resPEM3 <- rrisk.BayesPEM(x=x, n=n,k=k,
     prior.pi=c(1,1),prior.se=c(se.a,se.b),prior.sp=c(sp.a,sp.b),
     misclass="compare")
resPEM3@results

#------------------------------------------
# Example of PEM model (k=1)
#------------------------------------------
# informative priors -> convergence is o.k.
resPEM4 <- rrisk.BayesPEM(x=2,n=10,k=1,prior.se=c(12,22),
 prior.sp=c(22,55),prior.pi=c(1,1))
resPEM4@results

# non-informative priors -> convergence of 'pi' is not o.k.
resPEM5 <- rrisk.BayesPEM(x=2,n=10,k=1,prior.se=c(1,1),
 prior.sp=c(1,1),prior.pi=c(1,1))
resPEM5@results

# informative priors -> convergence is o.k., without invoking graphical
# diagnostic interface
resPEM6 <- rrisk.BayesPEM(x=2,n=10,k=1,prior.se=c(12,22),
 prior.sp=c(22,55),prior.pi=c(1,1))
resPEM6@results

#------------------------------------------
# Example of ZIP model
#------------------------------------------
# generate ZIP data
pi <- 0.01
n <- 200
lambda <- 3.5
zip.data <- rep(0,n)
zip.data[sample(1:n,n*pi,replace=FALSE)] <- rpois(n*pi,lambda=lambda)

# estimate using Bayes model for zero inflated data
resZIP1 <- rrisk.BayesZIP(data=zip.data, prior.lambda=c(0,100),prior.pi=c(1,1),
 burn=100,update=1000)
resZIP1@results

# estimate using Bayes model for zero inflated data without invoking
# graphical diagnostic interface
rrisk.BayesZIP(data=zip.data, prior.lambda=c(0,100),prior.pi=c(1,1),
 burn=100,update=1000,simulation=TRUE)

# compare with naive results ignoring ZIP model
pi.crude <- sum(zip.data>0)/n
lambda.crude <- mean(zip.data)
print(pi.crude)
print(lambda.crude)
resZIP1@results

#------------------------------------------
# Examples of GUI functions
#------------------------------------------
data <- rpois(30, 4)
res <- ZIPGUI(data)

mod <- PEMGUI()

#-----------------------------------------------------
# Creating an instance of the 'bayesmodelClass'
#-----------------------------------------------------
new("bayesmodelClass")
# end of donttest
 
````
