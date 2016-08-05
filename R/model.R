

################################################################################
################################################################################
#' @description test Bayesian PEM models provide the posterior distribution for the true prevalence (\code{pi}),
#' diagnostic sensitivity (\code{se}) and specificity (\code{sp}) for a given empirical prevalence
#' estimate using physically pooled samples (if \code{k>1}) and priors for the model parameters.
#' The misclassification parameters (\code{se} and \code{sp}) can be specified
#' at the level of the pool or individual level of testing. On the other side,
#' the function estimates the true prevalence based on the results
#' (\code{x/n}) of an application study with individual samples (if \code{k=1}) using a diagnostic test, for
#' which some prior information on sensitivity and specificity is available.
#'
#' @details The Bayesian model for estimation prevalence, sensitivity and specificity has
#' in BRugs/Winbugs syntax following form for misclassification at the pool-level
#' (\code{k>1} and \code{misclass="pool"})
#' \preformatted{model{
#'
#'        pi ~ dbeta(prior.pi[1],prior.pi[2])
#'
#'        se ~ dbeta(prior.se[1],prior.se[2])
#'
#'        sp ~ dbeta(prior.sp[1],prior.sp[2])
#'
#'        p.neg <- pow(1-pi,k)
#'
#'        p <- (1-p.neg)*se + p.neg*(1-sp)
#'
#'        x ~ dbin(p,n)
#'
#'      }}
#' for misclassifications at the individual level (\code{k>1} and \code{misclass="individual"})
#' \preformatted{model{
#'
#'        pi ~ dbeta(prior.pi[1],prior.pi[2])
#'
#'        se ~ dbeta(prior.se[1],prior.se[2])
#'
#'        sp ~ dbeta(prior.sp[1],prior.sp[2])
#'
#'        ap <- pi*se + (1-pi)*(1-sp)
#'
#'        p <- 1- pow(1-ap,k)
#'
#'        x ~ dbin(p,n)
#'
#'      }}
#' and for comparison (\code{k>1})
#' \preformatted{model{
#'
#'        pi1 ~ dbeta(prior.pi[1],prior.pi[2])
#'
#'        pi2 ~ dbeta(prior.pi[1],prior.pi[2])
#'
#'        se ~ dbeta(prior.se[1],prior.se[2])
#'
#'        sp ~ dbeta(prior.sp[1],prior.sp[2])
#'
#'        x1 <- x
#'
#'        x2 <- x
#'
#'        p.neg <- pow(1-pi1,k)
#'
#'        p.pos <- (1-p.neg)*se + p.neg*(1-sp)
#'
#'        x1 ~ dbin(p.pos,n)
#'
#'        ap <- pi2*se + (1-pi2)*(1-sp)
#'
#'        p <- 1- pow(1-ap,k)
#'
#'        x2 ~ dbin(p,n)
#'
#'      }}
#' The application data (\code{k=1}) has one degree of freedom while the underlying model
#' has three unknown parameters. Thus, the model is not identifiable and informative
#' priors on at least two model parameters are required. The Bayesian model for estimation
#' prevalence, sensitivity and specificity takes a form
#' \preformatted{model{
#'
#'      x ~ dbin(p,n)
#'
#'      p <- pi * se + (1-pi) * (1-sp)
#'
#'      se ~ dbeta(prior.se[1],prior.se[2])
#'
#'      sp ~ dbeta(prior.sp[1],prior.sp[2])
#'
#'      pi  ~ dbeta(prior.pi[1],prior.pi[2])
#'
#'    }}
#'
#' @name rrisk.BayesPEM
#' @aliases rrisk.BayesPEM
#' @title Bayesian Prevalence estimation under misclassification (PEM)
#' @usage rrisk.BayesPEM(x, n, k, simulation=FALSE, prior.pi, prior.se, prior.sp,
#'  misclass="pool",chains=3, burn=1000, update=10000,
#'  workdir=getwd(), plots=FALSE)
#' @param x scalar value for number of pools (\code{k>1}) or individual outcomes (\code{k=1}) with positive test result
#' @param n scalar value for number of pools tested (\code{k>1}) or the sample size in application study (\code{k=1})
#' @param k scalar value for number of individual samples physically combined into one pool;
#' set \code{k>1} for pooled sampling and \code{k=1} for individual sampling
#' @param simulation logical, value \code{TRUE} means the function will be called within any simulation routine,
#' in this case the graphical diagnostic interface will not be invoked (default \code{FALSE})
#' @param prior.pi numeric vector containing parameters of a beta distribution as prior for prevalence \code{pi}, e.g. \code{pi} \eqn{\sim} \code{prior.pi(*,*)=beta(*,*)}
#' @param prior.se numeric vector containing parameters of a beta distribution as prior for sensitivity \code{se}, e.g. \code{se} \eqn{\sim} \code{prior.se(*,*)=beta(*,*)}
#' @param prior.sp numeric vector containing parameters of a beta distribution as prior for specificity \code{sp}, e.g. \code{sp} \eqn{\sim} \code{prior.sp(*,*)=beta(*,*)}
#' @param misclass character with legal character entries \code{pool}, \code{individual} or \code{compare}; ignored if k=1
#' @param chains positive single numeric value, number of independent MCMC chains (default 3)
#' @param burn positive single numeric value, length of the burn-in period (default 1000)
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
#' \item{\code{update}}{length of update iterations for estimation}
#' @note The convergence of the model is assessed by the user using diagnostic plots
#' provided by the \pkg{BRugs} package.
# @seealso nothing...
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
                           prior.pi,
                           prior.se,
                           prior.sp,
                           misclass = "pool",
                           chains = 3,
                           burn = 1000,
                           update = 10000,
                           workdir = getwd(),
                           plots = FALSE
)
{

  checkInput(x,n,k, prior.pi, prior.se, prior.sp, chains, burn, update, misclass, workdir, plots)


  #-----------------------------------------------------------------------------
  # create a temporary directory
  #-----------------------------------------------------------------------------
  setwd(workdir)

    #-----------------------------------------------------------------------------
  # replace the out list below with S4 return objects ...
  #-----------------------------------------------------------------------------
  out<-new("bayesmodelClass")

  #-----------------------------------------------------------------------------
  # create nodes (parameters to be estimated) for the model and table for parameters of prior distribution for each node
  #-----------------------------------------------------------------------------

  # nodes <- c("pi","se","sp")
  # nodes.prior <- matrix(1,nrow=3, ncol=2)
  # colnames(nodes.prior) <- c("Parameter1","Parameter2")
  # rownames(nodes.prior) <- nodes
  # nodes.prior["pi",] <- prior.pi
  # nodes.prior["se",] <- prior.se
  # nodes.prior["sp",] <- prior.sp


  #-----------------------------------------------------------------------------
  # write model
  #-----------------------------------------------------------------------------
  if(k==1){
  } else if (k>1){
    if(misclass == "individual") {
      model_string <- function(pi_prior = c(1, 1), se_prior, sp_prior) {
        sprintf("model {
                  pi ~  dbeta(%g, %g)                 # prevalence
                  se ~  dbeta(%g, %g)                 # sensitivity
                  sp ~  dbeta(%g, %g)                 # specificity
                  ap <- pi*se + (1-pi)*(1-sp)         # apparent prevalence
                  x  ~  dbin(ap, n)                   # number of positive samples

                  #inits# pi, se, sp
                  #monitor# pi, se, sp
        }", pi_prior[1], pi_prior[2], se_prior[1], se_prior[2], sp_prior[1], sp_prior[2])
      }
    }else if(misclass == "individual-fix-sp") {
      model_string <- function(pi_prior, se_prior, sp_fix) {
        sprintf("model {
                  pi ~  dbeta(%g, %g)                 # prevalence
                  se ~  dbeta(%g, %g)                 # sensitivity
                  sp <- %g                            # fixed specificity
                  ap <- pi*se + (1-pi)*(1-sp)         # apparent prevalence
                  x  ~  dbin(ap, n)                   # number of positive samples

                  #inits# pi, se
                  #monitor# pi, se
        }", pi_prior[1], pi_prior[2], se_prior[1], se_prior[2], sp_fix)
      }
    }else if(misclass == "individual-fix-se") {
      model_string <- function(pi_prior, se_fix, sp_prior) {
        sprintf("model {
                pi ~  dbeta(%g, %g)                   # prevalence
                se ~  %g                              # sensitivity
                sp <- dbeta(%g, %g)                   # fixed specificity
                ap <- pi*se + (1-pi)*(1-sp)           # apparent prevalence
                x  ~  dbin(ap, n)                     # number of positive samples

                #inits# pi, se
                #monitor# pi, se
      }", pi_prior[1], pi_prior[2], se_fix, sp_prior[1], sp_prior[2])
      }
    } else if(misclass == "pool"){
      model_string <- function(pi_prior, seP_prior, spP_prior) {
        sprintf("model {
                pi  ~  dbeta(%g, %g)                 # prevalence
                seP ~  dbeta(%g, %g)                 # (pooled test) sensitivity
                spP ~  dbeta(%g, %g)                 # (pooled test) specificity
                p.neg <- pow(1-pi, k)                # probability of a disease-free pool
                p  <- (1-p.neg)*seP + p.neg*(1-spP)  # probability for a pool
                # to test positive
                x   ~  dbin(p, n)                    # number of positive pools

                #inits# pi, seP, spP
                #monitor# pi, seP, spP
      }", pi_prior[1], pi_prior[2], seP_prior[1], seP_prior[2], spP_prior[1], spP_prior[2])
      }
    }
  }








#-----------------------------------------------------------------------------
# output
#-----------------------------------------------------------------------------
#write.table(out$model,file="doc.txt",quote=FALSE,row.names=FALSE,col.names=FALSE)
return(out)
} # end of function


