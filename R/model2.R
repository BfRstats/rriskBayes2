
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
#'  misclass="pool",chains=3, burn=1000, thin=1, update=10000,
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
#' @param thin positive single numeric value (default 1). The samples from every kth iteration will be used for
#'        inference, where k is the value of thin. Setting \code{thin > 1} can help to reduce the autocorrelation
#'        in the sample.
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
#' @examples
#' \donttest{
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
#'      prior.pi=c(1,1),prior.se=c(se.a,se.b),prior.sp=c(sp.a,sp.b),
#'      misclass="individual")
#' resPEM1@@results
#'
#' # Estimate using Bayes model at pool level
#' resPEM2 <- rrisk.BayesPEM(x=x, n=n,k=k,
#'      prior.pi=c(1,1),prior.se=c(se.a,se.b),prior.sp=c(sp.a,sp.b),
#'      misclass="pool")
#' resPEM2@@results
#'
#' # Estimate using Bayes model compared
#' resPEM3 <- rrisk.BayesPEM(x=x, n=n,k=k,
#'      prior.pi=c(1,1),prior.se=c(se.a,se.b),prior.sp=c(sp.a,sp.b),
#'      misclass="compare")
#' resPEM3@@results
#'
#' #------------------------------------------
#' # Example of PEM model (k=1)
#' #------------------------------------------
#' # informative priors -> convergence is o.k.
#' resPEM4<-rrisk.BayesPEM(x=2,n=10,k=1,prior.se=c(12,22),
#'  prior.sp=c(22,55),prior.pi=c(1,1))
#' resPEM4@@results
#'
#' # non-informative priors -> convergence is not o.k.
#' resPEM5<-rrisk.BayesPEM(x=2,n=10,k=1,prior.se=c(1,1),
#'  prior.sp=c(1,1),prior.pi=c(1,1))
#' resPEM5@@results
#'
#' # informative priors -> convergence is o.k., without invoking
#' # graphical diagnostic interface
#' resPEM6<-rrisk.BayesPEM(x=2,n=10,k=1,prior.se=c(12,22),
#'  prior.sp=c(22,55),prior.pi=c(1,1))
#' resPEM6@@results
#' }

rrisk.BayesPEM <- function(x, n, k, simulation=FALSE,
    prior.pi, prior.se, prior.sp, misclass="pool",
    chains=3, burn=1000, thin=1, update=10000, workdir=getwd(), plots=FALSE)
{
  #-----------------------------------------------------------------------------
  # checking input
  #-----------------------------------------------------------------------------
  if(missing(x) | missing(n) | missing(k))
  { on.exit(return(invisible(NA)))
    stop("INVALID INPUT, missing one or more of the following arguments: 'x', 'n', 'k'!",call.=FALSE)
  }
  if(missing(prior.pi) | missing(prior.se) | missing(prior.sp))
  { on.exit(return(invisible(NA)))
    stop("INVALID INPUT, missing one or more of the following prior arguments: 'prior.pi', 'prior.se', 'prior.sp'!",call.=FALSE)
  }
  if(!is.numeric(chains) | !is.numeric(burn) | !is.numeric(update) | !is.numeric(x) | !is.numeric(n) | !is.numeric(k))
  { on.exit(return(invisible(NA)))
    stop("INVALID INPUT, one or more of the following arguments is not numeric: 'chains', 'burn', 'update'!",call.=FALSE)
  }
  max.l <- max(length(x),length(n),length(k)) # variable pool size disabled
  if((length(x) == max.l)*(length(n) == max.l)*(length(k) == max.l) != 1)
  { on.exit(return(invisible(NA)))
    stop("INVALID INPUT, equal length of arguments 'x', 'n' and 'k' required.",call.=FALSE)
  }
  if(chains<=0 | burn<=0 | update<=0 | thin<1)
  { on.exit(return(invisible(NA)))
    stop("INVALID INPUT, one or more of the following arguments is not positive: 'chains', 'burn', 'update', 'thin'!",call.=FALSE)
  }
  if(!is.element(misclass,c("pool","individual","compare")))
  { on.exit(return(invisible(NA)))
    stop("INVALID INPUT, misclass should be 'pool','individual' or 'compare'!",call.=FALSE)
  }
  try.setwd<-try(setwd(workdir),silent=TRUE)
  if(inherits(try.setwd, "try-error"))
  { on.exit(return(invisible(NA)))
    error.mess<-paste("INVALID INPUT, the working directory could not be found! Your input is \n",workdir)
    stop(error.mess,call.=FALSE)
  }
  if(!is.logical(plots))
  { on.exit(return(invisible(NA)))
    stop("INVALID INPUT, the argument 'plots' should be of type logical!",call.=FALSE)
  }
  #-----------------------------------------------------------------------------
  # plausibility check on priors
  #-----------------------------------------------------------------------------
  if(min(prior.pi)<=0 | min(prior.se)<=0 | min(prior.sp)<=0)
  { on.exit(return(invisible(NA)))
    stop("Parameters of the beta prior distribution for prevalence, sensitivity and specificity should be strictly positive!",call.=FALSE)
  }
  if(length(prior.pi)!=2 | length(prior.se)!=2 | length(prior.sp)!=2)
  { on.exit(return(invisible(NA)))
    stop("Two parameters for the beta prior distribution for prevalence, sensitivity and specificity are required",call.=FALSE)
  }

  #-----------------------------------------------------------------------------
  # create a temporary directory
  #-----------------------------------------------------------------------------
  setwd(workdir)
  #if(!file.exists("tempBayes"))
  #{
  #  suppressWarnings(dir.create("tempBayes"))
  #} else if(file.exists("tempBayes"))
  #{
  #  # system("rmdir tempBayes")
  #  unlink("tempBayes", recursive=TRUE)
  #  suppressWarnings(dir.create("tempBayes"))
  #}
  # setwd(file.path(workdir,"tempBayes")) # no more `cd`

  #-----------------------------------------------------------------------------
  # replace the out list below with S4 return objects ...
  #-----------------------------------------------------------------------------
  out<-new("bayesmodelClass")

  #-----------------------------------------------------------------------------
  # create nodes (parameters to be estimated) for the model and table for parameters of prior distribution for each node
  #-----------------------------------------------------------------------------

  nodes <- c("pi","se","sp")
  nodes.prior <- matrix(1,nrow=3, ncol=2)
  colnames(nodes.prior) <- c("Parameter1","Parameter2")
  rownames(nodes.prior) <- nodes
  nodes.prior["pi",] <- prior.pi
  nodes.prior["se",] <- prior.se
  nodes.prior["sp",] <- prior.sp

  #-----------------------------------------------------------------------------
  # write model
  #-----------------------------------------------------------------------------
  if(k==1)
  { cat(
    "model {
      se ~ dbeta(prior.se[1],prior.se[2])
      sp ~ dbeta(prior.sp[1],prior.sp[2])
      pi  ~ dbeta(prior.pi[1],prior.pi[2])
      x ~ dbin(p,n)
      p <- pi*se+(1-pi)*(1-sp)
    }",file="model.txt")
  } else if (k>1)
  { if(misclass == "pool")
    {  cat(
      "model{
        pi ~ dbeta(prior.pi[1],prior.pi[2])
        se ~ dbeta(prior.se[1],prior.se[2])
        sp ~ dbeta(prior.sp[1],prior.sp[2])
        p.neg <- pow(1-pi,k)
        p <- (1-p.neg)*se + p.neg*(1-sp)
        x ~ dbin(p,n)
      }",file="model.txt")
    } else if(misclass == "individual")
    { cat(
      "model{
        pi ~ dbeta(prior.pi[1],prior.pi[2])
        se ~ dbeta(prior.se[1],prior.se[2])
        sp ~ dbeta(prior.sp[1],prior.sp[2])
        ap <- pi*se + (1-pi)*(1-sp)
        p <- 1- pow(1-ap,k)
        x ~ dbin(p,n)
      }",file="model.txt")
    } else if(misclass == "compare")
    { cat(
      "model{
        pi1 ~ dbeta(prior.pi[1],prior.pi[2])
        pi2 ~ dbeta(prior.pi[1],prior.pi[2])
        se ~ dbeta(prior.se[1],prior.se[2])
        sp ~ dbeta(prior.sp[1],prior.sp[2])

        p.neg <- pow(1-pi1,k)
        p.pos <- (1-p.neg)*se + p.neg*(1-sp)
        x1 ~ dbin(p.pos,n)

        ap <- pi2*se + (1-pi2)*(1-sp)
        p <- 1- pow(1-ap,k)
        x2 ~ dbin(p,n)

        d <- pi1 - pi2

      }",file="model.txt")
     nodes <- c("pi1","pi2","se","sp","d")
    }
  }

  #-----------------------------------------------------------------------------
  # write data
  #-----------------------------------------------------------------------------
  if(k>1)
  {data.set <- list(x=x,n=n,k=k,prior.se=prior.se,prior.sp=prior.sp,prior.pi=prior.pi)
  } else if(k==1)
  { data.set <- list(x=x,n=n,prior.se=prior.se,prior.sp=prior.sp,prior.pi=prior.pi)
  }

  dput(data.set, file = "data.txt",control="keepNA")

  #-----------------------------------------------------------------------------
  # check model, load data and compile three chains
  #-----------------------------------------------------------------------------
  cat("--------------------------------------
      -----------------------------\n")
  cat("Begin model fitting...\n")
  #modelCheck("model.txt")
  jmod <- jags.model("model.txt", data=data.set, n.chains=chains)

  #modelData("data.txt")
  #read.bugsdata("data.txt")


  #gibt es nicht mehr in rjags
  #modelCompile(numChains = chains)
  #modelGenInits()

  #-----------------------------------------------------------------------------
  # run model with burn-in
  #-----------------------------------------------------------------------------
  # Burn-in:
  cat( "Burning in the MCMC chain...\n" )

  #modelUpdate(burn,thin=thin)
  #samplesSet(nodes)
  #modelUpdate(update,thin=thin)


  nIter <- ceiling((update*thin)/chains)
  update(obj=jmod, n.iter=burn)
  # The saved MCMC chain:
  cat( "Sampling final MCMC chain...\n" )
  codaSamples = coda.samples(jmod, variable.names=nodes, n.iter=nIter, thin=thin)
  # resulting codaSamples object has these indices:
  #   codaSamples[[ chainIdx ]][ stepIdx , paramIdx ]
  mcmcChain <- as.matrix(codaSamples)
  samplesSample <- as.data.frame(mcmcChain)


  # run model with iterations to reach stationary joint distribution

  #-----------------------------------------------------------------------------
  # ?berpr?fe, ob model erfolgreich gefittet wurde
  #-----------------------------------------------------------------------------
  plotDiag.check<-function(nodes){
    X11()
    plot1<-plotDiag(nodes,plotnumber=1)
    X11()
    plot2<-plotDiag(nodes,plotnumber=2)
    plot3<-plotDiag(nodes,plotnumber=3)
    output<-list(plot1=plot1,plot2=plot2,plot3=plot3)
  }
  try.result<-try(plotDiag.check(nodes),silent=FALSE)
  if (inherits(try.result, "try-error")) {
    if(!is.null(dev.list())) {graphics.off()}
    on.exit(return("ERROR"))
    stop("Error occured during model fitting...",call.=FALSE)
  }

  #-----------------------------------------------------------------------------
  # check diagnostic plots for convergence
  #-----------------------------------------------------------------------------
  if(!simulation)
  { out@convergence <- FALSE
    convergence <- diagnostics(nodes,plots)
    if(!is.logical(convergence)){ # falls "ERROR" zur?ckgeliefert
      on.exit(return(convergence))
      stop("Any error occured during model fitting...ABORD function",call.=FALSE)
    } else { # falls TRUE oder FALSE zur?ckgeliefert
      if(!convergence)
      {
        on.exit(return(invisible(NA)))
        # tidy-up
        # setwd(workdir)
        #unlink("tempBayes",recursive=TRUE)
        cat("End model fitting...\n")
        cat("-------------------------------------------------------------------\n")
        stop("Process has been cancelled by the user due to non-convergence",call.=FALSE)
      }
      out@convergence <- convergence
    }
  }

  #-----------------------------------------------------------------------------
  # estimation results from BRugs package
  #-----------------------------------------------------------------------------
  results <- .samplesStats1(nodes,beg = samplesGetBeg(), end = samplesGetEnd(),
    firstChain = samplesGetFirstChain(),
    lastChain = samplesGetLastChain())
  out@results <- results


  #-----------------------------------------------------------------------------
  # collect posterior joint distribution
  #-----------------------------------------------------------------------------
  if(misclass != "compare") pi.out <- samplesSample["pi"]
  if(misclass == "compare") {
    pi.out <- samplesSample["pi1"]
    pi.out <- samplesSample["pi2"]
    pi.out <- samplesSample["d"]
    }
  se.out <- samplesSample["se"]
  sp.out <- samplesSample["sp"]
  if(misclass == "compare") {
    out@jointpost <- data.frame(pi=NA,se=NA,sp=NA)
  } else {
    out@jointpost <- data.frame(pi=pi.out,se=se.out,sp=sp.out)
  }


  #-----------------------------------------------------------------------------
  # full model description (which can be run in Winbugs if necessary)
  #-----------------------------------------------------------------------------

  out@model <- paste("
  # pi=prevalence on element level
  # se=sensitivity of diagnostic method
  # sp=specificity of diagnostic method

  # model
  ",paste(suppressWarnings(read.table(file="model.txt",sep="\n",stringsAsFactors=FALSE))[,1],collapse="\n"),"
  # data
  ",paste(suppressWarnings(read.table(file="data.txt",sep="\n",stringsAsFactors=FALSE))[,1],collapse="\n"),"\n",collapse=" ")

  #-----------------------------------------------------------------------------
  # nodes names
  #-----------------------------------------------------------------------------
  out@chains <- chains
  out@nodes <- nodes
  out@burn<-burn
  out@update<-update
  cat("End model fitting...\n")
  cat("-------------------------------------------------------------------\n")

  #-----------------------------------------------------------------------------
  # tidy-up
  #-----------------------------------------------------------------------------
  file.remove("data.txt")
  file.remove("model.txt")
  # setwd(workdir)
  # unlink("tempBayes",recursive=TRUE)

  #-----------------------------------------------------------------------------
  # output
  #-----------------------------------------------------------------------------
  #write.table(out$model,file="doc.txt",quote=FALSE,row.names=FALSE,col.names=FALSE)
  return(out)
} # end of function




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
#'
#'    lambda ~ dunif(prior.lambda[1],prior.lambda[2])
#'
#'    pi ~ dbeta(prior.pi[1],prior.pi[2])
#'
#'    for (i in 1:n) {
#'
#'                 y[i]  ~ dpois(mu[i])
#'
#'                 mu[i] <- I[i] * lambda
#'
#'                I[i] ~ dbern(pi)
#'
#'    }
#'
#'  }}
#'
#' @name rrisk.BayesZIP
#' @aliases rrisk.BayesZIP
#' @title Bayes estimation of a zero inflated Poisson (ZIP) model
#' @usage rrisk.BayesZIP(data, prior.lambda=c(1,10), prior.pi=c(0.8,1), simulation=FALSE,
#'  chains=3, burn=1000, thin=1, update=10000, workdir=getwd(), plots=FALSE)
#' @param data matrix, data frame or data set with positive integers, including zeros and of the minimal length 10
#' @param prior.lambda numeric vector containing minimum and maximum of a uniform
#' distribution used as prior for the Poisson parameter \code{lambda}, e.g. \cr \code{lambda} \eqn{\sim} \code{prior.lambda(*,*)=unif(*,*)}
#' @param prior.pi numeric vector containing parameters of a beta distribution
#' describing prior knowledge about prevalence (proportion of contaminated samples), e.g. \cr \code{pi} \eqn{\sim} \code{prior.pi(*,*)=beta(*,*)}
#' @param simulation logical, value \code{TRUE} means the function will be called within any simulation routine,
#' in this case the graphical diagnostic interface will not be invoked (default \code{FALSE})
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
#'  burn=100,update=1000)
#' resZIP@@results
#'
#' # estimate using Bayes model for zero inflated data without invoking
#' # graphical diagnostic interface
#' rrisk.BayesZIP(data=zip.data, prior.lambda=c(0,100),prior.pi=c(1,1),
#'  burn=100,update=1000,simulation=TRUE)
#'
#' # compare with naive results ignoring ZIP model
#' pi.crude <- sum(zip.data>0)/n
#' lambda.crude <- mean(zip.data)
#' print(pi.crude)
#' print(lambda.crude)
#' resZIP@@results
#' }

rrisk.BayesZIP <- function(data, prior.lambda=c(1,10), prior.pi=c(0.8,1), simulation=FALSE,
  chains=3, burn=1000, thin=1, update=10000, workdir=getwd(), plots=FALSE)
{
  #-----------------------------------------------------------------------------
  # checking input
  #-----------------------------------------------------------------------------
  if(missing(data))
  { on.exit(return(invisible(NA)))
    stop("INVALID INPUT, missing one or more of the function arguments: 'data', 'prior.lambda', 'prior.pi'!",call.=FALSE)
  }
  if(!is.numeric(data)| !is.numeric(prior.lambda) | !is.numeric(prior.pi) | !is.numeric(chains) | !is.numeric(burn) | !is.numeric(update))
  { on.exit(return(invisible(NA)))
    stop("INVALID INPUT, one or more of the following arguments is not numeric: 'data', 'prior.lambda', 'prior.pi', 'chains', 'burn', 'update'!",call.=FALSE)
  }
  if(!is.logical(plots))
  { on.exit(return(invisible(NA)))
    stop("INVALID INPUT, the argument 'plots' should be of type logical!",call.=FALSE)
  }
  if(chains<=0 | burn<=0 | update<=0 | thin<1)
  { on.exit(return(invisible(NA)))
    stop("INVALID INPUT, one or more of the following arguments is not positive: 'chains', 'burn', 'update', 'thin'!",call.=FALSE)
  }
  if(length(prior.lambda)!=2 |length(prior.pi)!=2 )
  { on.exit(return(invisible(NA)))
    stop("INVALID INPUT, the arguments 'prior.lambda' and 'prior.pi' should be of length 2!",call.=FALSE)
  }
  try.setwd<-try(setwd(workdir),silent=TRUE)
  if(inherits(try.setwd, "try-error"))
  { on.exit(return(invisible(NA)))
    error.mess<-paste("INVALID INPUT, the working directory could not be found! Your input is \n",workdir)
    stop(error.mess,call.=FALSE)
  }
  #-----------------------------------------------------------------------------
  # plausibility check on data
  #-----------------------------------------------------------------------------
  if(length(data)<10)
  { on.exit(return(invisible(NA)))
    stop("Data set too small for this purpose!",call.=FALSE)
  }
  if(min(data)<0)
  { on.exit(return(invisible(NA)))
    stop("Negative counts in data set are not allowed!",call.=FALSE)
  }
  if(any(abs(round(data)-data)!=0))
  { on.exit(return(invisible(NA)))
    stop("Data set contains non-integer values!",call.=FALSE)
  }
  #-----------------------------------------------------------------------------
  # plausibility check on priors
  #-----------------------------------------------------------------------------
  if(min(prior.pi)<=0)
  { on.exit(return(invisible(NA)))
    stop("Parameters of the beta prior distribution for prevalence should not be negative!",call.=FALSE)
  }
  if(min(prior.lambda)<0)
  { on.exit(return(invisible(NA)))
    stop("Parameters of the uniform prior distribution used as prior for the Poisson parameter should be strictly positive!",call.=FALSE)
  }

  #-----------------------------------------------------------------------------
  # create a temporary directory
  #-----------------------------------------------------------------------------
  setwd(workdir)
#  if(!file.exists("tempBayes"))
#  {
#    suppressWarnings(dir.create("tempBayes"))
#  } else if(file.exists("tempBayes"))
#  {
#    # system("rmdir tempBayes")
#    unlink("tempBayes", recursive=TRUE)
#    suppressWarnings(dir.create("tempBayes"))
#  }
  # setwd(file.path(workdir,"tempBayes")) # no more `cd`

  #-----------------------------------------------------------------------------
  # replace the out list below with S4 return objects ...
  #-----------------------------------------------------------------------------
  out<-new("bayesmodelClass")
  #out <- list()

  #-----------------------------------------------------------------------------
  # write model
  #-----------------------------------------------------------------------------
  cat(
  "model{
    lambda ~ dunif(prior.lambda[1],prior.lambda[2])
    pi ~ dbeta(prior.pi[1],prior.pi[2])
    for (i in 1:n) {
    y[i]  ~ dpois(mu[i])
    mu[i] <- I[i] * lambda
    I[i] ~ dbern(pi)
    }
  }",file="model.txt")

  #-----------------------------------------------------------------------------
  # write data
  #-----------------------------------------------------------------------------
  data.set <- list(y=data, n=length(data), prior.lambda=prior.lambda,prior.pi=prior.pi)
  dput(data.set, file = "data.txt",control="keepNA")

  #-----------------------------------------------------------------------------
  # check model, load data and compile three chains
  #-----------------------------------------------------------------------------
  cat("-------------------------------------------------------------------\n")
  cat("Begin model fitting...\n")
  modelCheck("model.txt")
  modelData("data.txt")
  modelCompile(numChains = chains)

  #-----------------------------------------------------------------------------
  # generate inits
  #-----------------------------------------------------------------------------
  modelGenInits()

  #-----------------------------------------------------------------------------
  # run model with burn-in
  #-----------------------------------------------------------------------------
  modelUpdate(burn,thin=thin)
  samplesSet(c("pi","lambda"))
  # run model with iterations to reach stationary joint distribution
  modelUpdate(update,thin=thin)

  #-----------------------------------------------------------------------------
  # ?berpr?fe, ob model erfolgreich gefittet wurde
  #-----------------------------------------------------------------------------
  plotDiag.check<-function(nodes){
    X11()
    plot1<-plotDiag(nodes,plotnumber=1)
    X11()
    plot2<-plotDiag(nodes,plotnumber=2)
    plot3<-plotDiag(nodes,plotnumber=3)
    output<-list(plot1=plot1,plot2=plot2,plot3=plot3)
  }
  try.result<-try(plotDiag.check(nodes=c("pi","lambda")),silent=FALSE)
  if (inherits(try.result, "try-error")) {
    if(!is.null(dev.list())) {graphics.off()}
    on.exit(return("ERROR"))
    stop("Error occured during model fitting...",call.=FALSE)
  }

  #-----------------------------------------------------------------------------
  # check diagnostic plots for convergence
  #-----------------------------------------------------------------------------
  if(!simulation)
  { out@convergence <- FALSE
    convergence <- diagnostics(nodes=c("pi","lambda"),plots)
    if(!is.logical(convergence)){ # falls "ERROR" zur?ckgeliefert
      on.exit(return(convergence))
      stop("Any error occured during model fitting...ABORD function",call.=FALSE)
    } else { # falls TRUE oder FALSE zur?ckgeliefert
      if(!convergence)
      {
        on.exit(return(invisible(NA)))
        # tidy-up
        # setwd(workdir)
        #unlink("tempBayes",recursive=TRUE)
        cat("End model fitting...\n")
        cat("-------------------------------------------------------------------\n")
        stop("Process has been cancelled by the user due to non-convergence",call.=FALSE)
      }
      out@convergence <- convergence
    }
  }

  #-----------------------------------------------------------------------------
  # estimation results from BRugs package
  #-----------------------------------------------------------------------------
  results <- .samplesStats1(node=c("pi","lambda"),beg = samplesGetBeg(), end = samplesGetEnd(),
    firstChain = samplesGetFirstChain(),
    lastChain = samplesGetLastChain())
  out@results <- results

  #-----------------------------------------------------------------------------
  # collect posterior joint distribution
  #-----------------------------------------------------------------------------
  pi.out <- samplesSample["pi"]
  lambda.out <- samplesSample["lambda"]
  out@jointpost <- data.frame(pi=pi.out,lambda=lambda.out)

  #-----------------------------------------------------------------------------
  # full model description (which can be run in Winbugs if necessary)
  #-----------------------------------------------------------------------------
  out@model <- paste(" # lambda=Poisson density parameter
  # pi=prevalence of contaminated sampes
  #
  # model
  ",paste(suppressWarnings(read.table(file="model.txt",sep="\n",stringsAsFactors=FALSE))[,1],collapse="\n"),"
  # data
  ",paste(suppressWarnings(read.table(file="data.txt",sep="\n",stringsAsFactors=FALSE))[,1],collapse="\n"),"\n",collapse=" ")

  #-----------------------------------------------------------------------------
  # nodes names
  #-----------------------------------------------------------------------------
  out@nodes <- c("lambda","pi")
  out@burn<-burn
  out@update<-update
  out@chains<-chains
  cat("End model fitting...\n")
  cat("-------------------------------------------------------------------\n")

  #-----------------------------------------------------------------------------
  # tidy-up
  #-----------------------------------------------------------------------------
  file.remove("data.txt")
  file.remove("model.txt")
  #  setwd(workdir)
  #  unlink("tempBayes",recursive=TRUE)

  #-----------------------------------------------------------------------------
  # output
  #-----------------------------------------------------------------------------
  #write.table(out$model,file="doc.txt",quote=FALSE,row.names=FALSE,col.names=FALSE)
  return(out)
} # end of function rrisk.BayesZIP()






