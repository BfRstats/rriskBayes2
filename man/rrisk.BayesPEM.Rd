% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/models.R
\name{rrisk.BayesPEM}
\alias{rrisk.BayesPEM}
\title{Bayesian Prevalence Estimation under Misclassification (PEM)}
\usage{
rrisk.BayesPEM(x, n, k, simulation = FALSE, 
 prior.pi = c(1, 1), prior.se, prior.sp,
 misclass = "pool", chains = 3, burn = 4000, 
 thin = 1, update = 5000, 
 max.time = "15minutes", plots = FALSE)
}
\arguments{
\item{x}{Scalar value for number of pools (\code{k>1}) or individual outcomes 
(\code{k=1}) with positive test result.}

\item{n}{Scalar value for number of pools tested (\code{k>1}) or the sample 
size in application study (\code{k=1}).}

\item{k}{Scalar value for number of individual samples physically combined 
into one pool;
set \code{k>1} for pooled sampling and \code{k=1} for individual sampling (default 1).}

\item{prior.pi}{Numeric vector containing parameters of a beta distribution 
as prior for prevalence \code{pi}, e.g., \cr \code{pi ~ prior.pi(*,*) = beta(*,*)}.}

\item{prior.se}{Numeric vector containing parameters of a beta distribution 
as prior for sensitivity \code{se}, e.g., \cr \code{se ~ prior.se(*,*) = beta(*,*)}. 
For fixed sensitivity scalar value.}

\item{prior.sp}{Numeric vector containing parameters of a beta distribution 
as prior for specificity \code{sp}, e.g., \cr \code{sp ~ prior.sp(*,*) = beta(*,*)}. 
For fixed specifity scalar value.}

\item{simulation}{Not used any longer.}

\item{misclass}{Character with legal character entries: \cr
\code{individual}, \code{individual-fix-se}, \code{individual-fix-sp}, \code{individual-fix-se-sp}, \cr
\code{pool}, \code{pool-fix-se}, \code{pool-fix-sp} or \code{pool-fix-se-sp}.\cr
\code{fix-se}: fixed sensitivity\cr
\code{fix-sp}: fixed specifity\cr
\code{fix-se-sp}: fixed sensitivity AND fixed specifity.}

\item{chains}{Positive single numeric value, number of independent MCMC 
chains (default 3).}

\item{burn}{Positive single numeric value, length of the burn-in period 
(default 4000).}

\item{thin}{Positive single numeric value (default 1). The samples from every 
k-th iteration will be used for inference, where k is the value of thin. 
Setting thin > 1 can help to reduce the autocorrelation in the sample.}

\item{update}{Positive single numeric value, length of update iterations for estimation (default 5000).}

\item{max.time}{Maximum time for which the function is allowed to extend the chains. Acceptable units include 'seconds', 'minutes', 'hours', 'days', 'weeks' (default "15minutes") (see \link[runjags]{autorun.jags}).}

\item{plots}{Logical, if \code{TRUE} the diagnostic plots will be displayed.}
}
\value{
The function \code{rrisk.BayesZIP} returns an instance of the \code{\linkS4class{bayesmodelClass}}
class containing following information:
\item{\code{convergence}}{Logical, whether the model has converged (assessed by the user).}
\item{\code{results}}{Data frame containing statitsics of the posterior distribution.}
\item{\code{jointpost}}{Data frame giving the joint posterior probability distribution.}
\item{\code{nodes}}{Names of the parameters jointly estimated by the Bayes model.}
\item{\code{model}}{Model in rjags/JAGS syntax as a character string.}
\item{\code{chains}}{Number of independent MCMC chains.}
\item{\code{burn}}{Length of burn-in period.}
\item{\code{update}}{Length of update iterations for estimation.}
}
\description{
Bayesian PEM models provide the posterior distribution for the 
true prevalence (\code{pi}), diagnostic sensitivity (\code{se}) and 
specificity (\code{sp}) for a given empirical prevalence estimate using 
physically pooled samples (if \code{k>1}) and priors for the model parameters.
The misclassification parameters (\code{se} and \code{sp}) can be specified
at the level of the pool or individual level of testing. On the other side,
the function estimates the true prevalence based on the results
(\code{x/n}) of an application study with individual samples (if \code{k=1}) 
using a diagnostic test, for
which some prior information on sensitivity and/or specificity is available.
}
\details{
The application data (\code{k=1}) has one degree of freedom while the 
underlying model has three unknown parameters. Thus, the model is not 
identifiable and informative priors on at least two model parameters are required. The Bayesian models for estimation prevalence, sensitivity and specificity take the following forms in rjags/JAGS (originally BRugs/Winbugs) syntax:

 \cr
Misclassifications at the individual level (\code{k=1} and \code{misclass="individual"})
\preformatted{model{

      pi  ~ dbeta(prior.pi[1], prior.pi[2])

      se ~ dbeta(prior.se[1], prior.se[2])

      sp ~ dbeta(prior.sp[1], prior.sp[2])

      ap <- pi * se + (1-pi) * (1-sp)
      
      x ~ dbin(ap, n)
   }}
Misclassifications at the individual level with fixed sensitivity 
resp. specifity (\code{k=1} and \code{misclass="individual-fix-sp"})
Both se and sp could be set to fixed values.
\preformatted{model{
      pi ~ dbeta(prior.pi[1], prior.pi[2])

      se resp. sp ~ dbeta(prior.sp[1], prior.sp[2])

      sp resp. se <-  fix

      ap <- pi*se + (1-pi)*(1-sp)

      x ~ dbin(p,n)
   }}
For misclassification at the pool-level (\code{k>1} and \code{misclass="pool"})
\preformatted{model{

      pi ~ dbeta(prior.pi[1], prior.pi[2])

      se ~ dbeta(prior.se[1], prior.se[2])

      sp ~ dbeta(prior.sp[1], prior.sp[2])

      p.neg <- pow(1-pi, k)

      p <- (1-p.neg)*se + p.neg*(1-sp)

      x ~ dbin(p, n)
     }}
For misclassification at the pool-level (\code{k>1}) and fixed 
     sensitivity and specifity (\code{misclass="pool-fix-se-sp"})
\preformatted{model{

      pi ~ dbeta(prior.pi[1], prior.pi[2])

      se <- fix

      sp <- fix

      p.neg <- pow(1-pi,k)

      p <- (1-p.neg)*se + p.neg*(1-sp)

      x ~ dbin(p,n)
     }}
and for comparison (\code{k>1})
\preformatted{model{

      pi1 ~ dbeta(prior.pi[1], prior.pi[2])

      pi2 ~ dbeta(prior.pi[1], prior.pi[2])

      se ~ dbeta(prior.se[1], prior.se[2])

      sp ~ dbeta(prior.sp[1], prior.sp[2])

      x1 <- x

      x2 <- x

      p.neg <- pow(1-pi1, k)

      p.pos <- (1-p.neg)*se + p.neg*(1-sp)

      x1 ~ dbin(p.pos, n)

      ap <- pi2*se + (1-pi2)*(1-sp)

      p <- 1- pow(1-ap,k)

      x2 ~ dbin(p,n)
     }}
}
\examples{
\donttest{
------------------------------------------
Example of PEM model
------------------------------------------
# generate PEM data at individual level

n <- 100
x <- 14
k <- 1
pi_prior <- c(1, 1)
se_prior <- c(64, 4)
sp_prior <- c(94, 16)
misclass <- "individual"

# run model

resPEM <- rrisk.BayesPEM(x = x, n = n, k = k, 
prior.pi = pi_prior, prior.se = se_prior, prior.sp = sp_prior,  
misclass = misclass, plots=TRUE)
  
# generate data for pooled sampling and fixed sensitivity and specifity

n <- 100
x <- 14
k <- 4
pi_prior <- c(1, 1)
sp_prior <- 1
se_prior <- 0.7
misclass <- "pool-fix-se-sp"

# run model

resPEM <- rrisk.BayesPEM(x = x, n = n, k = k, 
prior.pi = pi_prior, prior.se = se_prior, prior.sp = sp_prior, 
misclass = misclass, plots=TRUE)
}
}
\references{
Greiner, M., Belgoroski, N., NN (2011). Estimating prevalence 
using diagnostic data with pooled samples: R function \code{rrisk.BayesPEM}. 
J.Stat.Software (in preparation).

Cowling, D.W., Gardner, I.A. and Johnson W.O. (1999). Comparison 
of methods for estimation of individual-level prevalence based on pooled 
samples, Prev.Vet.Med. 39: 211-225.
\cr
\cr
Rogan, W.J. and B. Gladen (1978). Estimating prevalence from the results of 
a screening test. Am. J. Epidemiol. 107: 71-76.
}
\keyword{manip}

