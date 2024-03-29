
################################################################################
#' This packages provides a collection of functions for fitting Bayesian models
#' (related to the \code{rrisk} project). The functions can be used as stand-alone
#' applications or launched during an \code{rrisk} session.
#' \cr
#' The following Bayesian models are implemented in this package:
#' \describe{
#' \item{PEM}{prevalence estimation under misclassification}
#' \item{ZIP}{estimation of a zero inflated Poisson model}
#' \item{ZINB}{zero inflated negative binomial model}}
#' 
#' This package is a part of the \code{rrisk} project and contains functions for
#' fitting Bayesian models using the R package \pkg{rjags} and the software \pkg{JAGS}. This package does not
#' depend on the whole \code{rrisk} project and can be used separately. The \code{rrisk}
#' project can be downloaded from \url{http://www.bfr.bund.de/cd/52158}.
#'
#' @name rriskBayes-package
#' @aliases rriskBayes
#' @docType package
#' @concept Bayes model, Rogan-Gladen, prevalence, sensitivity, specificity, MCMC,
#' ZIP, Poisson, zero-inflation, rrisk, stat-up
#' @title Predefined Bayes models fitted with Markov chain Monte-Carlo (MCMC)
#' @author Natalia Belgorodski \email{belgorodski@@stat-up.de} (\acronym{STAT-UP} Statistical Consulting, Germany), \cr
#' Matthias Greiner \email{matthias.greiner@@bfr.bund.de} (Federal Institute for Risk Assessment, Germany), \cr
#' Alexander Engelhardt \email{engelhardt@@stat-up.de} (\acronym{STAT-UP} Statistical Consulting, Germany), \cr
#' Christina Yassouridis \email{christina-yassouridis@@stat-up.com} (\acronym{STAT-UP} Statistical Consulting, Germany)
#' @note See also the documentation to the \acronym{R}-package \pkg{\link{rjags}}.
#' @keywords package
#' @examples
#' \donttest{
#' 
#' library(rriskBayes2)
#' 
#' #------------------------------------------
#' # Example of PEM model (k>1)
#' #------------------------------------------
#' pi <- 0.01
#' se <- 0.96
#' se.n <- 1000
#' sp <- 0.99
#' sp.n <- 1000
#' n <- sample(10:1000,1,replace=TRUE)  # stochatsic sample size
#' k <- sample(5:50,1,replace=FALSE)    # stochastic pool size
#' 
#' # Parameters for beta priors
#' se.a <- se.n*se+1
#' se.b <- se.n*(1-se)+1
#' sp.a <- sp.n*sp+1
#' sp.b <- sp.n*(1-sp)+1
#' 
#' # Random number of positive pools (x) considering uncertainty of se and sp
#' ap <- pi*se + (1-pi)*(1-sp)
#' p.pos <- 1-(1-ap)^k
#' x <- rbinom(1,prob=p.pos,size=n)
#' 
#' # Estimate using Bayes model at individual level
#' resPEM1 <- rrisk.BayesPEM(x=x, n=n,k=k,
#'                           prior.pi=c(1,1),prior.se=c(se.a,se.b),prior.sp=c(sp.a,sp.b),
#'                           misclass="pool")
#' resPEM1@results
#' 
#' # Estimate using Bayes model at pool level
#' resPEM2 <- rrisk.BayesPEM(x=x, n=n, k=k,
#' prior.pi=c(1,1),prior.se=c(se.a,se.b),prior.sp=c(sp.a,sp.b),
#' misclass="pool")
#' resPEM2@results
#' 
#' # Estimate using Bayes model compared
#' n <- sample(10:100, 1)
#' k <- sample(2:10, 1)
#' x <- sample(1:n, 1)
#' pi_prior <- c(1, 1)
#' se_prior <- c(64, 4)
#' sp_prior <- c(94, 16)
#' misclass <- "compare"
#' 
#' resPEM3 <- rrisk.BayesPEM(x=x, n=n,k=k,
#'                           prior.pi=c(1,1),prior.se=se_prior, prior.sp=sp_prior,
#'                           misclass="compare")
#' resPEM3@results
#' 
#' #------------------------------------------
#' # Example of PEM model (k=1)
#' #------------------------------------------
#' # informative priors -> convergence is o.k.
#' resPEM4 <- rrisk.BayesPEM(x=2,n=10,k=1,prior.se=c(12,22),
#' prior.sp=c(22,55),prior.pi=c(1,1))
#' resPEM4@results
#' 
#' # non-informative priors -> convergence of 'pi' is not o.k.
#' resPEM5 <- rrisk.BayesPEM(x=2,n=10,k=1,prior.se=c(1,1),
#'                           prior.sp=c(1,1),prior.pi=c(1,1))
#' resPEM5@results
#' 
#' # informative priors -> convergence is o.k., without invoking graphical
#' # diagnostic interface
#' resPEM6 <- rrisk.BayesPEM(x=2,n=10,k=1,prior.se=c(12,22),
#'                           prior.sp=c(22,55),prior.pi=c(1,1))
#' resPEM6@results
#' 
#' #------------------------------------------
#' # Example of ZIP model
#' #------------------------------------------
#' # generate ZIP data
#' pi <- 0.01
#' n <- 200
#' lambda <- 3.5
#' zip.data <- rep(0,n)
#' zip.data[sample(1:n,n*pi,replace=FALSE)] <- rpois(n*pi,lambda=lambda)
#' 
#' # estimate using Bayes model for zero inflated data
#' resZIP1 <- rrisk.BayesZIP(data=zip.data, prior.lambda=c(0,100),prior.pi=c(1,1),
#'                           burn=100,update=4000)
#' resZIP1@results
#' 
#' # estimate using Bayes model for zero inflated data without invoking
#' # graphical diagnostic interface
#' rrisk.BayesZIP(data=zip.data, prior.lambda=c(0,100),prior.pi=c(1,1),
#'                burn=100,update=4000,simulation=TRUE)
#' 
#' # compare with naive results ignoring ZIP model
#' pi.crude <- sum(zip.data>0)/n
#' lambda.crude <- mean(zip.data)
#' print(pi.crude)
#' print(lambda.crude)
#' resZIP1@results
#' #------------------------------------------
#' # Examples of GUI functions
#' #------------------------------------------
#' data <- rpois(30, 4)
#' res <- ZIPGUI(data)
#' 
#' mod <- PEMGUI()
#' 
#' #------------------------------------------
#' # Example of ZINB model
#' #------------------------------------------
#' #generate ZINB data
#' pi <- 0.1
#' n <- 200
#' zinb.data <- rep(0,n)
#' zinb.data[sample(1:n, n*pi, replace = FALSE)] <- rpois(n*pi, lambda = 4)
#' 
#' # estimate using Bayes model for zero inflated data
#' resZINB <- rrisk.BayesZINB(data = zinb.data, prior.pi = c(1,1),
#'                            burn = 100, update = 4000, max.time = '40seconds')
#' resZINB
#' 
#' #-----------------------------------------------------
#' # Creating an instance of the 'bayesmodelClass'
#' #-----------------------------------------------------
#' new("bayesmodelClass")
#' # end of donttest
#' 
#' }

NULL