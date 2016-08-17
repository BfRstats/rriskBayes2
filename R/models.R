################################################################################
################################################################################
#' @description Bayesian PEM models provide the posterior distribution for the true prevalence (\code{pi}),
#' diagnostic sensitivity (\code{se}) and specificity (\code{sp}) for a given empirical prevalence
#' estimate using physically pooled samples (if \code{k>1}) and priors for the model parameters.
#' The misclassification parameters (\code{se} and \code{sp}) can be specified
#' at the level of the pool or individual level of testing. On the other side,
#' the function estimates the true prevalence based on the results
#' (\code{x/n}) of an application study with individual samples (if \code{k=1}) using a diagnostic test, for
#' which some prior information on sensitivity and specificity is available.
#'
#' @details 
#' The Bayesian model for estimation prevalence, sensitivity and specificity has
#' in BRugs/Winbugs syntax following form.
#' The application data (\code{k=1}) has one degree of freedom while the underlying model
#' has three unknown parameters. Thus, the model is not identifiable and informative
#' priors on at least two model parameters are required. The Bayesian model for estimation
#' prevalence, sensitivity and specificity takes a form
#' \preformatted{model{
#'
#'       pi  ~ dbeta(prior.pi[1],prior.pi[2])
#'
#'       se ~ dbeta(prior.se[1],prior.se[2])
#'
#'       sp ~ dbeta(prior.sp[1],prior.sp[2])
#'
#'       ap <- pi * se + (1-pi) * (1-sp)
#'       
#'       x ~ dbin(ap,n)
#'    }}
#' for misclassifications at the individual level (\code{k=1} and \code{misclass="individual"})
#' \preformatted{model{
#' 
#'       pi ~ dbeta(1, 1)
#'
#'       se ~ dbeta(prior.se[1],prior.se[2])
#'
#'       sp ~ dbeta(prior.sp[1],prior.sp[2])
#'
#'       ap <- pi*se + (1-pi)*(1-sp)
#'
#'       x ~ dbin(ap,n)
#'    }}
#' for misclassifications at the individual level with fixed sensitivity resp. specifity (\code{k=1} and \code{misclass="individual-fix-sp"})
#' \preformatted{model{
#'        pi ~ dbeta(prior.pi[1], prior.pi[2])
#'
#'        se resp. sp ~ dbeta(prior.sp[1], prior.sp[2])
#'
#'        sp resp. se ~ fix
#'
#'        ap <- pi*se + (1-pi)*(1-sp)
#'
#'        x ~ dbin(p,n)
#'    }}
#' For misclassification at the pool-level (\code{k>1} and \code{misclass="pool"})
#' \preformatted{model{
#'
#'       pi ~ dbeta(prior.pi[1], prior.pi[2])
#'
#'       se ~ dbeta(prior.se[1], prior.se[2])
#'
#'       sp ~ dbeta(prior.sp[1], prior.sp[2])
#'
#'       p.neg <- pow(1-pi,k)
#'
#'       p <- (1-p.neg)*se + p.neg*(1-sp)
#'
#'       x ~ dbin(p,n)
#'      }}
#' and for comparison (\code{k>1})
#' \preformatted{model{
#'
#'       pi1 ~ dbeta(prior.pi[1], prior.pi[2])
#'
#'       pi2 ~ dbeta(prior.pi[1], prior.pi[2])
#'
#'       se ~ dbeta(prior.se[1], prior.se[2])
#'
#'       sp ~ dbeta(prior.sp[1], prior.sp[2])
#'
#'       x1 <- x
#'
#'       x2 <- x
#'
#'       p.neg <- pow(1-pi1, k)
#'
#'       p.pos <- (1-p.neg)*se + p.neg*(1-sp)
#'
#'       x1 ~ dbin(p.pos, n)
#'
#'       ap <- pi2*se + (1-pi2)*(1-sp)
#'
#'       p <- 1- pow(1-ap,k)
#'
#'       x2 ~ dbin(p,n)
#'      }}
#' @name rrisk.BayesPEM
#' @aliases rrisk.BayesPEM
#' @title Bayesian Prevalence estimation under misclassification (PEM)
#' @usage rrisk.BayesPEM(x, n, k, simulation=FALSE, 
#'  prior.pi, prior.se, prior.sp,
#'  misclass="pool", chains=3, burn=1000, thin=1, update=10000,
#'  workdir=getwd(), plots=FALSE)
#' @param x scalar value for number of pools (\code{k>1}) or individual outcomes (\code{k=1}) with positive test result
#' @param n scalar value for number of pools tested (\code{k>1}) or the sample size in application study (\code{k=1})
#' @param k scalar value for number of individual samples physically combined into one pool;
#' set \code{k>1} for pooled sampling and \code{k=1} for individual sampling
#' @param simulation not used any longer
#' @param prior.pi numeric vector containing parameters of a beta distribution as prior for prevalence \code{pi}, e.g. \code{pi} ~ \code{prior.pi(*,*)=beta(*,*)}
#' @param prior.se numeric vector containing parameters of a beta distribution as prior for sensitivity \code{se}, e.g. \code{se} ~ \code{prior.se(*,*)=beta(*,*)}. For fixed sensitivity scalar value.
#' @param prior.sp numeric vector containing parameters of a beta distribution as prior for specificity \code{sp}, e.g. \code{sp} ~ \code{prior.sp(*,*)=beta(*,*)}. For fixed specifity scalara value.
#' @param misclass character with legal character entries \code{pool}, \code{individual} or \code{individual-fix-se} or \code{individual-fix-sp}.
#' @param chains positive single numeric value, number of independent MCMC chains (default 3)
#' @param burn positive single numeric value, length of the burn-in period (default 1000)
#' @param thin positive single numeric value (default 1). The samples from every kth iteration will be used for inference, where k is the value of thin. Setting thin > 1 can help to reduce the autocorrelation in the sample.
#' @param update positive single numeric value, length of update iterations for estimation (default 10000)
#' @param workdir character string giving working directory to store temporary data (default \code{getwd()})
#' @param plots logical, if \code{TRUE} the diagnostic plots will be displayed in separate windows
#' @return The function \code{rrisk.BayesPEM} returns an instance of the \code{\linkS4class{bayesmodelClass}}
#' class containing following informations
#' \item{\code{convergence}}{logical, whether the model has converged (assessed by the user)}
#' \item{\code{results}}{data frame containing statistics of the posterior distribution}
#' \item{\code{jointpost}}{data frame giving the joint posterior probability distribution}
#' \item{\code{nodes}}{names of the parameters jointly estimated by the Bayes model}
#' \item{\code{model}}{model in BRugs/Winbugs syntax as a character string}
#' \item{\code{chains}}{number of independent MCMC chains}
#' \item{\code{burn}}{length of burn-in period}
#' \item{\code{thin}}{the thinning interval to be used}
#' \item{\code{update}}{length of update iterations for estimation}
#' @keywords manip
#' @export
# @references Greiner, M., Belgoroski, N., NN (2011). Estimating prevalence using diagnostic data
# with pooled samples: R function \code{rrisk.BayesPEM}. J.Stat.Software (in preparation).
#' @references Cowling, D.W., I.A. Gardner and W.O. Johnson (1999). Comparison of methods for estimation
#'   of individual-level prevalence based on pooled samples, Prev.Vet.Med. 39: 211-225.
#' \cr
#' \cr
#' Rogan, W.J. and B. Gladen (1978). Estimating prevalence from the results of a screening test. Am. J. Epidemiol. 107: 71-76.
#' @examples TO WRITE


rrisk.BayesPEM <- function(x,
                           n,
                           k,
                           prior.pi = c(1, 1),
                           prior.se,
                           prior.sp,
                           misclass = "pool",
                           chains = 3,
                           burn = 1000,
                           thin = 1,
                           update = 10000,
                           workdir = getwd(),
                           plots = FALSE
)
{
  # -----------------------------------------------------------------------------
  # check input arguments
  # -----------------------------------------------------------------------------
  checkInputPEM(x, n, k, prior.pi, prior.se, prior.sp, chains, burn, thin, update, misclass, workdir, plots)

  #-----------------------------------------------------------------------------
  # create outlist
  #-----------------------------------------------------------------------------
  out<-new("bayesmodelClass")

  #-----------------------------------------------------------------------------
  # a priori model definitions
  #-----------------------------------------------------------------------------
  #data
  if (misclass == "pool")
    jags_data <- list(n = n, x = x, k = k)
  else if(misclass == "compare")
    jags_data <- list(n = n, x1 = x, x2 = x, k = k)
  else
    jags_data <- list(n = n, x = x)
  
  #wrapper function for inits_function with only one argument chain as required by autorun.jags
  inits <- function(chain) inits_functionPEM(chain, misclass)
  
  #define model
  model_function <- modelFunctionPEM(misclass)
  model <- model_function(prior.pi, prior.se, prior.sp)

  #-----------------------------------------------------------------------------
  # run model
  #-----------------------------------------------------------------------------
  jags_res <- autorun.jags(
    model    = model,
    data     = jags_data,
    n.chains = chains,
    inits    = inits,
    startburnin = burn,
    startsample = update,
    max.time = "3m",
    method   = "rjags",
    thin     = thin,
    plots    = FALSE
  )

  if(plots)
    plotDiag(jags_res)
  
  #-----------------------------------------------------------------------------
  # output
  #-----------------------------------------------------------------------------
  out@convergence <- checkPSRF(jags_res)
  out@nodes <- jags_res$monitor
  out@model <- writeModelPEM(misclass)
  out@chains <- chains
  out@burn <- burn
  out@update <- update
  out@jointpost <- (sample(jags_res))$mcmc[[1]]
  out@results <- jags_res

return(out)
} # end of function rrisk.BayesPEM





################################################################################
################################################################################
#' @description Zero-inflated Poisson data are count data with an excess number of zeros. The
#' ZIP model involves the Poisson parameter \code{lambda} and the prevalence
#' parameter \code{pi}.
#'
#' @details The ZIP model applies to count data and can be interpreted as a mixture
#' distribution with one component comprising the 'true' zeros and another component
#' of Poisson distributed values with density parameter \code{lambda}. The prevalence
#' parameter \code{pi} refers to the proportion of the second, true non-zero
#' component.
#' \cr \cr
#' The Bayesian model for estimation prevalence and lambda parameter has
#' in BRugs/Winbugs syntax following form 
#' \preformatted{model{
#'    lambda ~ dunif(prior.lambda[1], prior.lambda[2])
#'
#'    pi ~ dbeta(prior.pi[1], prior.pi[2]) 
#'
#'    for (i in 1:n) {  
#'                    y[i]  ~ dpois(mu[i])
#' 
#'                    mu[i] <- I[i] * lambda  
#'
#'                    I[i] ~ dbern(pi)  
#'                    }   
#'  }
#'}
#'
#' @name rrisk.BayesZIP
#' @aliases rrisk.BayesZIP
#' @title Bayes estimation of a zero inflated Poisson (ZIP) model
#' @usage rrisk.BayesZIP(data, prior.lambda=c(1,10), prior.pi=c(0.8,1), simulation=FALSE,
#'  chains=3, burn=1000, thin=1, update=10000, workdir=getwd(), plots=FALSE)
#' @param data matrix, data frame or data set with positive integers, including zeros and of the minimal length 10
#' @param prior.lambda numeric vector containing minimum and maximum of a uniform
#' distribution used as prior for the Poisson parameter \code{lambda}, e.g. \cr \code{lambda} ~ \code{prior.lambda(*,*)=unif(*,*)}
#' @param prior.pi numeric vector containing parameters of a beta distribution
#' describing prior knowledge about prevalence (proportion of contaminated samples), e.g. \cr \code{pi} ~ \code{prior.pi(*,*)=beta(*,*)}
#' @param simulation not used any longer
#' @param chains positive single numeric value, number of independent MCMC chains (default 3)
#' @param burn positive single numeric value, length of the burn-in period (default 1000)
#' @param thin positive single numeric value (default 1). The samples from every kth iteration will be used for 
#'        inference, where k is the value of thin. Setting \code{thin > 1} can help to reduce the autocorrelation
#'        in the sample.
#' @param update positive single numeric value, length of update iterations for estimation (default 10000)
#' @param workdir character string giving working directory to store temporary data (default \code{getwd()})
#' @param plots logical, if \code{TRUE} the diagnostic plots will be displayed in separate windows 
#' @return The function \code{rrisk.BayesZIP} returns an instance of the \code{\linkS4class{bayesmodelClass}}
#' class containing following informations
#' \item{\code{convergence}}{logical, whether the model has converged (assessed by the user)}
#' \item{\code{results}}{data frame containing statitsics of the posterior distribution}
#' \item{\code{jointpost}}{data frame giving the joint posterior probability distribution}
#' \item{\code{nodes}}{names of the parameters jointly estimated by the Bayes model}
#' \item{\code{model}}{model in BRugs/Winbugs syntax as a character string}
#' \item{\code{chains}}{number of independent MCMC chains}
#' \item{\code{burn}}{length of burn-in period}
#' \item{\code{update}}{length of update iterations for estimation}
#' @note The convergence of the model should be checked using the diagnostic plots
#' see the package \pkg{BRugs}, see also \pkg{zicounts}.
# @seealso nothing...
#' @keywords manip
#' @export
#' @references Bohning, D., E. Dietz, P. Schlattman, L. Mendonca, and U. Kirchner (1999). 
#' The zero-inflated Poisson model and the decayed, missing and filled teeth index in 
#' dental epidemiology. Journal of the Royal Statistical Society, Series A 162, 195-209.
#' @examples
#' \donttest{
#' #------------------------------------------
#' # Example of ZIP model
#' #------------------------------------------
#' # generate ZIP data
#' pi<-0.01
#' n<-200
#' lambda<-3.5
#' zip.data<-rep(0,n)
#' zip.data[sample(1:n,n*pi,replace=FALSE)]<-rpois(n*pi,lambda=lambda)
#'
#' # estimate using Bayes model for zero inflated data
#' resZIP<-rrisk.BayesZIP(data=zip.data, prior.lambda=c(0,100),prior.pi=c(1,1),
#'  burn=100,update=4000)
#' resZIP@@results
#'
#' # compare with naive results ignoring ZIP model
#' pi.crude <- sum(zip.data>0)/n
#' lambda.crude <- mean(zip.data)
#' print(pi.crude)
#' print(lambda.crude)
#' resZIP@@results
#' }

rrisk.BayesZIP <-
  function(data,
           prior.lambda = c(1, 10),
           prior.pi = c(0.8, 1),
           simulation = NULL,
           chains = 3,
           burn = 1000,
           thin = 1,
           update = 10000,
           workdir = getwd(),
           plots = FALSE)
  {
    # -----------------------------------------------------------------------------
    # check input arguments
    # -----------------------------------------------------------------------------
    checkInputZIP(data, prior.lambda, prior.pi, chains, burn, thin, update, workdir, plots)
    
    #-----------------------------------------------------------------------------
    # create outlist
    #-----------------------------------------------------------------------------
    out <- new("bayesmodelClass")
    
    #-----------------------------------------------------------------------------
    # a priori model definitions
    #-----------------------------------------------------------------------------
    #data
    jags_data <- list(y = data, n = length(data))
    
    #wrapper function for inits_function with only one argument chain as required by autorun.jags
    inits <- function(chains) inits_functionZIP(chain = chains, data = jags_data)
    
    #define model
    model_function <- modelFunctionZIP
    model <- model_function(pi_prior = prior.pi, lambda_prior = prior.lambda)
    
    #-----------------------------------------------------------------------------
    # run model
    #-----------------------------------------------------------------------------
    jags_res <- autorun.jags(
      model    = model,
      data     = jags_data,
      n.chains = chains,
      inits    = inits,
      startburnin = burn,
      startsample = update,
      max.time = "3m",
      method   = "rjags",
      thin     = thin,
      plots    = FALSE
    )
    
    if (plots)
      plotDiag(jags_res)
    
    #-----------------------------------------------------------------------------
    # output
    #-----------------------------------------------------------------------------
    out@convergence <- checkPSRF(jags_res)
    out@nodes <- jags_res$monitor
    out@model <- writeModelZIP()
    out@chains <- chains
    out@burn <- burn
    out@update <- update
    out@jointpost <- (sample(jags_res))$mcmc[[1]]
    out@results <- jags_res
    
    return(out)
  } # end of function rrisk.BayesZIP()


################################################################################
################################################################################


#' @name rrisk.BayesZINB
#' @aliases rrisk.BayesZINB
#' @title Bayes estimation of a zero inflated negative binomial (ZINB) model
#' 
#' @export
rrisk.BayesZINB <-
  function(data,
           prior.pi = c(0.8, 1),
           simulation = NULL,
           chains = 3,
           burn = 1000,
           thin = 1,
           update = 10000,
           workdir = getwd(),
           plots = FALSE)
  {
    
  # -----------------------------------------------------------------------------
  # check input arguments
  # -----------------------------------------------------------------------------
  checkInputZINB(data, prior.pi, simulation, chains, burn, thin, update, workdir, plots)
  
 #-----------------------------------------------------------------------------
  # create outlist
  #-----------------------------------------------------------------------------
  out <- new("bayesmodelClass")
  
  #-----------------------------------------------------------------------------
  # a priori model definitions
  #-----------------------------------------------------------------------------
  #data
  jags_data <- list(y = data, n = length(data))
  
  #wrapper function for inits_function with only one argument chain as required by autorun.jags
  inits <- function(chains) inits_functionZINB(chain = chains, data = jags_data)
  
  #define model
  model_function <- modelFunctionZINB
  model <- model_function(pi_prior = prior.pi)

  #-----------------------------------------------------------------------------
  # run model
  #-----------------------------------------------------------------------------
  jags_res <- autorun.jags(
    model    = model,
    data     = jags_data,
    n.chains = chains,
    inits    = inits,
    startburnin = burn,
    startsample = update,
    max.time = "3m",
    method   = "rjags",
    thin     = thin,
    plots    = FALSE
  )
  
  #-----------------------------------------------------------------------------
  # plots
  #-----------------------------------------------------------------------------
  if (plots)
    plotDiag(jags_res)
  
  #-----------------------------------------------------------------------------
  # output
  #-----------------------------------------------------------------------------
  out@convergence <- checkPSRF(jags_res)
  out@nodes <- jags_res$monitor
  out@model <- writeModelZIP()
  out@chains <- chains
  out@burn <- burn
  out@update <- update
  out@jointpost <- (sample(jags_res))$mcmc[[1]]
  out@results <- jags_res
  
 return(jags_res)
} # end of function rrisk.BayesZIP()



