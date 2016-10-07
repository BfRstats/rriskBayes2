
# rrisk.BayesPEM ----------------------------------------------------------
#' @description Bayesian PEM models provide the posterior distribution for the 
#' true prevalence (\code{pi}), diagnostic sensitivity (\code{se}) and 
#' specificity (\code{sp}) for a given empirical prevalence estimate using 
#' physically pooled samples (if \code{k>1}) and priors for the model parameters.
#' The misclassification parameters (\code{se} and \code{sp}) can be specified
#' at the level of the pool or individual level of testing. On the other side,
#' the function estimates the true prevalence based on the results
#' (\code{x/n}) of an application study with individual samples (if \code{k=1}) 
#' using a diagnostic test, for
#' which some prior information on sensitivity and/or specificity is available.
#' @details 
#' The Bayesian model for estimation prevalence, sensitivity and specificity has
#' in BRugs/Winbugs syntax the following form:
#' The application data (\code{k=1}) has one degree of freedom while the 
#' underlying model has three unknown parameters. Thus, the model is not 
#' identifiable and informative priors on at least two model parameters are 
#' required. The Bayesian model for estimation prevalence, sensitivity and 
#' specificity takes the following forms: \cr
#' Misclassifications at the individual level (\code{k=1} and \code{misclass="individual"})
#' \preformatted{model{
#'
#'       pi  ~ dbeta(prior.pi[1], prior.pi[2])
#'
#'       se ~ dbeta(prior.se[1], prior.se[2])
#'
#'       sp ~ dbeta(prior.sp[1], prior.sp[2])
#'
#'       ap <- pi * se + (1-pi) * (1-sp)
#'       
#'       x ~ dbin(ap, n)
#'    }}
#' Misclassifications at the individual level with fixed sensitivity 
#' resp. specifity (\code{k=1} and \code{misclass="individual-fix-sp"})
#' Both se and sp could be set to fixed values.
#' \preformatted{model{
#'        pi ~ dbeta(prior.pi[1], prior.pi[2])
#'
#'        se resp. sp ~ dbeta(prior.sp[1], prior.sp[2])
#'
#'        sp resp. se <-  fix
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
#'       p.neg <- pow(1-pi, k)
#'
#'       p <- (1-p.neg)*se + p.neg*(1-sp)
#'
#'       x ~ dbin(p, n)
#'      }}
#' For misclassification at the pool-level (\code{k>1}) and fixed 
#'      sensitivity and specifity (\code{misclass="pool-fix-se-sp"})
#' \preformatted{model{
#'
#'       pi ~ dbeta(prior.pi[1], prior.pi[2])
#'
#'       se <- fix
#'
#'       sp <- fix
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
#' @usage rrisk.BayesPEM(x, n, k, simulation = FALSE, 
#'  prior.pi = c(1, 1), prior.se, prior.sp,
#'  misclass = "pool", chains = 3, burn = 4000, 
#'  thin = 1, update = 5000, 
#'  max.time = '15minutes', plots = FALSE)
#' @param x scalar value for number of pools (\code{k>1}) or individual outcomes 
#' (\code{k=1}) with positive test result
#' @param n scalar value for number of pools tested (\code{k>1}) or the sample 
#' size in application study (\code{k=1})
#' @param k scalar value for number of individual samples physically combined 
#' into one pool;
#' set \code{k>1} for pooled sampling and \code{k=1} for individual sampling (default 1)
#' @param simulation not used any longer
#' @param prior.pi numeric vector containing parameters of a beta distribution 
#' as prior for prevalence \code{pi}, e.g. \cr \code{pi ~ prior.pi(*,*) = beta(*,*)}
#' @param prior.se numeric vector containing parameters of a beta distribution 
#' as prior for sensitivity \code{se}, e.g. \cr \code{se ~ prior.se(*,*) = beta(*,*)}. 
#' For fixed sensitivity scalar value.
#' @param prior.sp numeric vector containing parameters of a beta distribution 
#' as prior for specificity \code{sp}, e.g. \cr \code{sp ~ prior.sp(*,*) = beta(*,*)}. 
#' For fixed specifity scalar value.
#' @param misclass character with legal character entries \cr
#' \code{individual}, \code{individual-fix-se}, \code{individual-fix-sp}, \code{individual-fix-se-sp}, \cr
#' \code{pool}, \code{pool-fix-se}, \code{pool-fix-sp} or \code{pool-fix-se-sp}.\cr
#' \code{fix-se}: fixed sensitivity\cr
#' \code{fix-sp}: fixed specifity\cr
#' \code{fix-se-sp}: fixed sensitivity and fixed specifity
#' @param chains positive single numeric value, number of independent MCMC 
#' chains (default 3)
#' @param burn positive single numeric value, length of the burn-in period 
#' (default 4000)
#' @param thin positive single numeric value (default 1). The samples from every 
#' k-th iteration will be used for inference, where k is the value of thin. 
#' Setting thin > 1 can help to reduce the autocorrelation in the sample.
#' @param update positive single numeric value, length of update iterations for 
#' estimation (default 5000)
#' @param max.time the maximum time for which the function is allowed to extend the chains. Acceptable units include 'seconds', 'minutes', 'hours', 'days', 'weeks' (default '15minutes') (see \link[runjags]{autorun.jags})
#' @param plots logical, if \code{TRUE} the diagnostic plots will be displayed
#' @return The function \code{rrisk.BayesPEM} returns an instance of the 
#' \code{\linkS4class{bayesmodelClass}} class containing following informations
#' \item{\code{convergence}}{logical, whether the model has converged 
#' (assessed by the user)}
#' \item{\code{results}}{data frame containing statistics of the posterior 
#' distribution}
#' \item{\code{jointpost}}{data frame giving the joint posterior probability 
#' distribution}
#' \item{\code{nodes}}{names of the parameters jointly estimated by the Bayes 
#' model}
#' \item{\code{model}}{model in BRugs/Winbugs syntax as a character string}
#' \item{\code{chains}}{number of independent MCMC chains}
#' \item{\code{burn}}{length of burn-in period}
#' \item{\code{thin}}{the thinning interval to be used}
#' \item{\code{update}}{length of update iterations for estimation}
#' @keywords manip
#' @export
#' @references Greiner, M., Belgoroski, N., NN (2011). Estimating prevalence 
#' using diagnostic data with pooled samples: R function \code{rrisk.BayesPEM}. 
#' J.Stat.Software (in preparation).
#' @references Cowling, D.W., Gardner, I.A. and Johnson W.O. (1999). Comparison 
#' of methods for estimation of individual-level prevalence based on pooled 
#' samples, Prev.Vet.Med. 39: 211-225.
#' \cr
#' \cr
#' Rogan, W.J. and B. Gladen (1978). Estimating prevalence from the results of 
#' a screening test. Am. J. Epidemiol. 107: 71-76.
#' @examples
#' \donttest{
#'------------------------------------------
#' Example of PEM model
#'------------------------------------------
#' # generate PEM data at individual level
#' 
#' n <- 100
#'x <- 14
#' k <- 1
#' pi_prior <- c(1, 1)
#' se_prior <- c(64, 4)
#' sp_prior <- c(94, 16)
#' misclass <- "individual"
#' 
#' # run model
#' 
#'resPEM <- rrisk.BayesPEM(x = x, n = n, k = k, 
#'prior.pi = pi_prior, prior.se = se_prior, prior.sp = sp_prior,  
#'misclass = misclass, plots=TRUE)
#'   
#' # generate data for pooled sampling and fixed sensitivity and specifity
#' 
#' n <- 100
#' x <- 14
#' k <- 4
#' pi_prior <- c(1, 1)
#' sp_prior <- 1
#' se_prior <- 0.7
#' misclass <- "pool-fix-se-sp"
#' 
#' # run model
#' 
#' resPEM <- rrisk.BayesPEM(x = x, n = n, k = k, 
#' prior.pi = pi_prior, prior.se = se_prior, prior.sp = sp_prior, 
#' misclass = misclass, plots=TRUE)
#' }
rrisk.BayesPEM <- function(x, n, k = 1,
                           prior.pi = c(1, 1), prior.se, prior.sp,
                           simulation = FALSE,
                           misclass = "pool",
                           chains = 3, burn = 4000, thin = 1, update = 5000,
                           max.time = '15minutes',
                           plots = FALSE
                           )
{
  ##### check input arguments #####
  checkInputPEM(x, n, k, prior.pi, prior.se, prior.sp, chains, burn, thin, 
                update, misclass, plots)

  ##### create outlist #####
  out <- new("bayesmodelClass")

  ##### a priori model definitions #####
  #data
  if (misclass == "pool" | misclass == "pool-fix-se" | misclass == "pool-fix-sp" 
      | misclass == "pool-fix-se-sp")
    jags_data <- list(n = n, x = x, k = k)
  else if(misclass == "compare")
    jags_data <- list(n = n, x1 = x, x2 = x, k = k)
  else
    jags_data <- list(n = n, x = x)
  
  #####initialization #####
  #wrapper function for inits_function with only one argument chain as required 
  #by autorun.jags
  inits <- function(chain) inits_functionPEM(chain, misclass = misclass)

  #define model
  model_function <- modelDefinitionPEM(misclass)
  model <- model_function(prior.pi, prior.se, prior.sp)

  startburnin = burn - 1000 #burn - adapt
  
  ##### run model #####
  jags_res <- autorun.jags(
    model    = model,
    data     = jags_data,
    n.chains = chains,
    inits    = inits,
    startburnin = startburnin,
    startsample = update,
    max.time = max.time,
    method   = "rjags",
    thin     = thin,
    plots    = FALSE
  )
  
  ##### check diagnostic plots for convergence #####
  if(!simulation)
  { 
    convergence <- diagnostics(jags_res, plots)
 
      if(!convergence)
      {
        on.exit(return(invisible(NA)))
        cat("End model fitting...\n")
        cat("-----------------------------------------------------------------\n")
        stop("Process has been cancelled by the user due to non-convergence",
             call.=FALSE)
      }
      out@convergence <- convergence
  }

  ##### output #####
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


# rrisk.BayesZIP ----------------------------------------------------------

#' @description Zero-inflated Poisson data are count data with an excess number 
#' of zeros. The ZIP model involves the Poisson parameter \code{lambda} and the 
#' prevalence parameter \code{pi}.
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
#' 
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
#' @usage rrisk.BayesZIP(data, prior.lambda = c(1, 10), prior.pi = c(0.8, 1), 
#'  simulation = FALSE, chains = 3, burn = 4000, 
#'  thin = 1, update = 5000, 
#'  max.time = '15minutes', plots = FALSE)
#' @param data matrix, data frame or data set with positive integers, including 
#' zeros and of the minimal length 10
#' @param prior.lambda numeric vector containing minimum and maximum of a uniform
#' distribution used as prior for the Poisson parameter \code{lambda}, e.g.
#' \code{lambda} ~ \code{prior.lambda(*,*)=unif(*,*)}
#' @param prior.pi numeric vector containing parameters of a beta distribution
#' describing prior knowledge about prevalence (proportion of contaminated 
#' samples), e.g. \code{pi} ~ \code{prior.pi(*,*)=beta(*,*)}
#' @param simulation not used any longer 
#' @param chains positive single numeric value, number of independent MCMC chains 
#' (default 3)
#' @param burn positive single numeric value, length of the burn-in period 
#' (default 4000)
#' @param thin positive single numeric value (default 1). The samples from every 
#' k-th iteration will be used for inference, where k is the value of thin. 
#' Setting \code{thin > 1} can help to reduce the autocorrelation in the sample.
#' @param update positive single numeric value, length of update iterations for 
#' estimation (default 5000)
#' @param max.time the maximum time for which the function is allowed to extend the chains. Acceptable units include 'seconds', 'minutes', 'hours', 'days', 'weeks' (default '15minutes') (see \link[runjags]{autorun.jags})
#' @param plots logical, if \code{TRUE} the diagnostic plots will be displayed 
#' in separate windows 
#' @return The function \code{rrisk.BayesZIP} returns an instance of the 
#' \code{\linkS4class{bayesmodelClass}}
#' class containing following informations
#' \item{\code{convergence}}{logical, whether the model has converged (assessed 
#' by the user)}
#' \item{\code{results}}{data frame containing statitsics of the posterior 
#' distribution}
#' \item{\code{jointpost}}{data frame giving the joint posterior probability 
#' distribution}
#' \item{\code{nodes}}{names of the parameters jointly estimated by the Bayes model}
#' \item{\code{model}}{model in BRugs/Winbugs syntax as a character string}
#' \item{\code{chains}}{number of independent MCMC chains}
#' \item{\code{burn}}{length of burn-in period}
#' \item{\code{update}}{length of update iterations for estimation}
#' @note The convergence of the model should be checked using the diagnostic plots
#' @keywords manip
#' @export
#' @references Bohning, D., Dietz, E., Schlattman, P., Mendonca,  L. and  Kirchner, U. (1999). 
#' The zero-inflated Poisson model and the decayed, missing and filled teeth index in 
#' dental epidemiology. Journal of the Royal Statistical Society, Series A 162:195-209.
#' @examples
#' \donttest{
#' #------------------------------------------
#' # Example of ZIP model
#' #------------------------------------------
#' # generate ZIP data
#' pi <- 0.01
#' n <- 200
#' lambda <- 3.5
#' zip.data <- rep(0,n)
#' zip.data[sample(1:n,n*pi,replace=FALSE)]<-rpois(n*pi,lambda=lambda)
#'
#' # estimate using Bayes model for zero inflated data
#' resZIP <- rrisk.BayesZIP(data = zip.data, 
#' prior.lambda = c(0,100), 
#' prior.pi = c(1,1),
#' burn = 4000, 
#' max.time = '40seconds',
#' update = 4000)
#' resZIP
#'
#' # compare with naive results ignoring ZIP model
#' pi.crude <- sum(zip.data>0)/n
#' lambda.crude <- mean(zip.data)
#' print(pi.crude)
#' print(lambda.crude)
#' resZIP@results
#' }

rrisk.BayesZIP <-  function(data,
                            prior.lambda = c(1, 10), prior.pi = c(0.8, 1),
                            simulation = FALSE,
                            chains = 3,
                            burn = 4000,
                            thin = 1,
                            update = 5000,
                            max.time = '15minutes',
                            plots = FALSE)
  {
    #####check input arguments #####
    checkInputZIP(data, prior.lambda, prior.pi, chains, burn, thin, update, plots)
    
    ##### create outlist #####
    out <- new("bayesmodelClass")
    
    ##### a priori model definitions #####
    #data
    jags_data <- list(y = data, n = length(data))
    
    #wrapper function for inits_function with only one argument chain as required by autorun.jags
    inits <- function(chains) inits_functionZIP(chain = chains, data = jags_data)
    
    #define model
    model_function <- modelDefinitionZIP
    model <- model_function(pi_prior = prior.pi, lambda_prior = prior.lambda)
    
    startburnin = burn - 1000 #burn - adapt
    
    ##### run model #####
    jags_res <- autorun.jags(
      model    = model,
      data     = jags_data,
      n.chains = chains,
      inits    = inits,
      startburnin = startburnin,
      startsample = update,
      max.time = max.time,
      method   = "rjags",
      thin     = thin,
      plots    = FALSE
    )
    
    ##### check diagnostic plots for convergence #####
    
    if(!simulation)
    { 
      convergence <- diagnostics(jags_res, plots)
      
      if(!convergence)
      {
        on.exit(return(invisible(NA)))
        cat("End model fitting...\n")
        cat("-----------------------------------------------------------------\n")
        stop("Process has been cancelled by the user due to non-convergence",
             call.=FALSE)
      }
      out@convergence <- convergence
    }
    
    ##### output #####
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


# rrisk.BayesZINB ---------------------------------------------------------

#' @description Zero-inflated Negative-Binomial data are count data with an excess number of zeros. The
#' ZINB model involves the prevalence parameter \code{pi}. The negative binomial distribution can be seen
#' as a poisson(\eqn{\lambda}) distribution, where \eqn{\lambda} follows a gamma distribution.
#' @details The ZINB model applies to count data and can be interpreted as a mixture distribution with one component
#' comprising the 'true' zeros and another component of negative-binomially distributed values with density parameter 
#' \eqn{\lambda} following a gamma distribution with shape and mean parameters
#' modelled as \code{dgamma(shape = 0.01, lambda = 0.01)}. The prevalence parameter \code{pi} refers to the proportion of the second, true 
#' non-zero component.
#' \cr \cr
#' The Bayesian model for estimation prevalence and lambda parameter has
#' in BRugs/Winbugs syntax following form
#' \preformatted{model{
#' 
#'    pi  ~ dbeta(prior.pi[1], prior.pi[2])
#'     
#'    dam ~ dgamma(0.01, 0.01)
#'    
#'    db ~ dgamma(0.01, 0.01)
#'
#'    for (i in 1:n) {
#'            y[i] ~ dpois(mu[i])
#'                    
#'            mu[i] <- I[i] * lambda[i]
#'                    
#'            I[i] ~ dbern(pi)
#'                    
#'            lambda[i] ~ dgamma(dam, db)
#'                    }
#'  }
#'}
#' @name rrisk.BayesZINB
#' @aliases rrisk.BayesZINB
#' @title Bayes estimation of a zero inflated negative binomial (ZINB) (also referred to as 'gamma-poisson') model
#' @usage rrisk.BayesZINB(data, prior.pi = c(0.8, 1), 
#' simulation = FALSE, chains = 3, burn = 4000,
#'  thin = 1, update = 5000,
#'   max.time = '15minutes', plots = FALSE)
#' @param data matrix, data frame or data set with positive integers, including zeros and of the minimal length 10
#' @param prior.pi numeric vector containing parameters of a beta distribution
#' describing prior knowledge about prevalence (proportion of contaminated samples), e.g. \cr \code{pi} ~ \code{prior.pi(*,*)=beta(*,*)}
#' @param simulation not used any longer
#' @param chains positive single numeric value, number of independent MCMC chains (default 3)
#' @param burn positive single numeric value, length of the burn-in period (default 4000)
#' @param thin positive single numeric value (default 1). The samples from every k-th iteration will be used for
#'        inference, where k is the value of thin. Setting \code{thin > 1} can help to reduce the autocorrelation
#'        in the sample.
#' @param update positive single numeric value, length of update iterations for estimation (default 5000)
#' @param max.time the maximum time for which the function is allowed to extend the chains. Acceptable units include 'seconds', 'minutes', 'hours', 'days', 'weeks' (default '15minutes') (see \link[runjags]{autorun.jags})
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
#' @export
#' @references Lunn, D. et al. (2012). The BUGS book: A practical introduction to Bayesian analysis. CRC press. p.353, 117
#' @examples
#' \donttest{
#' #------------------------------------------
#' #Example of a ZINB model (compare with rrisk.BayesZIP)
#' #------------------------------------------
#' #generate ZINB data
#' pi <- 0.1
#' n <- 200
#' zinb.data <- rep(0,n)
#' zinb.data[sample(1:n, n*pi, replace = FALSE)] <- rpois(n*pi, lambda = 4)
#'
#' # estimate using Bayes model for zero inflated data
#' resZINB <- rrisk.BayesZINB(data = zinb.data, prior.pi = c(1,1),
#'  burn = 100, update = 4000, max.time = '40seconds')
#' resZINB
#' 
#' #------------------------------------------
#' #Example of a ZINB model 
#' #------------------------------------------
#'set.seed(42)

#'n_true_neg <- 60
#'n_true_pos <- 33
#'n <- n_true_pos + n_true_neg
#'
#'prev_true <- n_true_pos / n
#'
#'a <- 6
#'b <- 2
#'lambda_true <- rgamma(n_true_pos, a, b)
#'
#'y_neg <- rep(0, n_true_neg)
#'y_pos <- rpois(n_true_pos, lambda_true)
#'y <- c(y_pos, y_neg)
#'
#'pi_prior     <- c(1, 1)
#'
#'resZINB <- rrisk.BayesZINB(data = y,
#'                             prior.pi = pi_prior,
#'                             simulation = FALSE,
#'                             chains = 3,
#'                             burn = 4000,
#'                             thin = 1,
#'                             max.time = '60seconds',
#'                             update = 10000,
#'                             plots = TRUE
#'                             )
#' }

rrisk.BayesZINB <-  function(data,
                             prior.pi = c(0.8, 1),
                             simulation = FALSE,
                             chains = 3,
                             burn = 4000,
                             thin = 1,
                             update = 5000,
                             max.time = '15minutes',
                             plots = FALSE)
  {
    
 ##### check input arguments #####
 checkInputZINB(data, prior.pi, chains, burn, thin, update, plots)
  
 ##### create outlist #####
  out <- new("bayesmodelClass")
  
  ##### a priori model definitions #####
  #data
  jags_data <- list(y = data, n = length(data))
  
  #wrapper function for inits_function with only one argument chain as required by autorun.jags
  inits <- function(chains) inits_functionZINB(chain = chains, data = jags_data)
  
  #define model
  model_function <- modelDefinitionZINB
  model <- model_function(pi_prior = prior.pi)

  startburnin = burn - 1000 #burn - adapt
  
  ##### run model #####
  jags_res <- autorun.jags(
    model    = model,
    data     = jags_data,
    n.chains = chains,
    inits    = inits,
    startburnin = startburnin,
    startsample = update,
    max.time = max.time,
    method   = "rjags",
    thin     = thin,
    plots    = FALSE
  )
  
  ##### check diagnostic plots for convergence #####
  if(!simulation)
  { 
    convergence <- diagnostics(jags_res, plots)
    
    if(!convergence)
    {
      on.exit(return(invisible(NA)))
      cat("End model fitting...\n")
      cat("-----------------------------------------------------------------\n")
      stop("Process has been cancelled by the user due to non-convergence",
           call.=FALSE)
    }
    out@convergence <- convergence
  }
  
  ##### output #####
  out@convergence <- checkPSRF(jags_res)
  out@nodes <- jags_res$monitor
  out@model <- writeModelZINB()
  out@chains <- chains
  out@burn <- burn
  out@update <- update
  out@jointpost <- (sample(jags_res))$mcmc[[1]]
  out@results <- jags_res
  
 return(out)
} # end of function rrisk.BayesZINB()



