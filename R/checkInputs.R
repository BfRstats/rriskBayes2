################################################################################
################################################################################

#-----------------------------------------------------------------------------
# common input for all models
#-----------------------------------------------------------------------------

checkMCparams <- function(chains, burn, thin, update){
  if( chains > 6)
  { on.exit(return(invisible(NA)))
    stop("The maximum number of chains is 5!",call.=FALSE)
  }
  
  if(chains<=0 | burn<=0 | update<=0 | thin <1)
  { on.exit(return(invisible(NA)))
    stop("INVALID INPUT, one or more of the following arguments is not positive: 
         'chains', 'burn', 'update', 'thin'!",call.=FALSE)
  }
  
  if(!is.numeric(chains) | !is.numeric(burn) | !is.numeric(update))
  { on.exit(return(invisible(NA)))
    stop("INVALID INPUT, one or more of the following arguments is not numeric: 
         'chains', 'burn', 'update', 'x', 'n' !",call.=FALSE)
  }
}

checkEnvironParams <- function(workdir, plots){
  
  try.setwd<-try(setwd(workdir),silent=TRUE)
  if(inherits(try.setwd, "try-error"))
  { on.exit(return(invisible(NA)))
    error.mess<-paste("INVALID INPUT, the working directory could not be found!
                      Your input is \n",workdir)
    stop(error.mess,call.=FALSE)
  }
  
  if(!is.logical(plots))
  { on.exit(return(invisible(NA)))
    stop("INVALID INPUT, the argument 'plots' should be of type logical!",
         call.=FALSE)
  }
}

################################################################################
################################################################################
#-----------------------------------------------------------------------------
# input for specific models
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
# input for rrisk.BayesPEM
#-----------------------------------------------------------------------------

checkInputPEM <- function(x, n, k,
                          prior.pi, prior.se, prior.sp,
                          chains, burn, thin, update,
                          misclass,
                          workdir, plots
                          )
{
  checkMCparams(chains, burn, thin, update)
  checkEnvironParams(workdir, plots)
  
  
  if (missing(x) | missing(n))
  { on.exit(return(invisible(NA)))
    stop("INVALID INPUT, missing one or more of the following arguments: 'x', 
         'n'!", call. = FALSE)
  }

  max.l <- max(length(x),length(n),length(k)) # variable pool size disabled
  if((length(x) == max.l)*(length(n) == max.l)*(length(k) == max.l) != 1)
  { on.exit(return(invisible(NA)))
    stop("INVALID INPUT, equal length of arguments 'x', 'n' and 'k' required.",
         call.=FALSE)
  }
  
  #-----------------------------------------------------------------------------
  # plausibility check on priors
  #-----------------------------------------------------------------------------

  if(missing(prior.pi) | missing(prior.se) | missing(prior.sp))
  { on.exit(return(invisible(NA)))
    stop("INVALID INPUT, missing one or more of the following prior arguments: 
         'prior.pi', 'prior.se', 'prior.sp'!", call. = FALSE)
  }
  
  if(min(prior.pi)<=0 | min(prior.se)<=0 | min(prior.sp)<=0)
  { on.exit(return(invisible(NA)))
    stop("Parameters of the beta prior distribution for prevalence, sensitivity
         and specificity should be strictly positive!",call.=FALSE)
  }
  
  #-----------------------------------------------------------------------------
  # plausibility check on method
  #-----------------------------------------------------------------------------
  
  if(!is.element(misclass,c("individual","individual-fix-sp", "individual-fix-se",
                            "individual-fix-se-sp", "pool", "pool-fix-se", 
                            "pool-fix-sp", "pool-fix-se-sp", "compare")))
  { on.exit(return(invisible(NA)))
    stop("INVALID INPUT, misclass should be 'individual', 'individual-fix-sp', 
         'individual-fix-se', 'indivdual-fix-se-sp', 'pool', 'pool-fix-se', 
         'pool-fix-sp', 'pool-fix-se-sp',  or 'compare'!",call.=FALSE)
  }
  
  if (misclass == "pool") {
    if(missing(k) | !is.numeric(k))
    { on.exit(return(invisible(NA)))
      stop("INVALID INPUT, missing or non-numeric 'k'!", call. = FALSE)
    }
    
  } else if (misclass == "individual-fix-sp" | misclass == "pool-fix-sp") {
    if (length(prior.sp) != 1) 
    { on.exit(return(invisible(NA)))
      stop("Specifity must be a fix value!")
    }
    
  } else if (misclass == "individual-fix-se" | misclass == "pool-fix-se") {
    if (length(prior.se) != 1) 
    { on.exit(return(invisible(NA)))
      stop("Sensitivity must be a fix value!")
    }
  } else if (misclass == "individual-fix-se-sp" | misclass == "pool-fix-se-sp") {
    if (length(prior.se) != 1 | length(prior.sp) != 1) 
    { on.exit(return(invisible(NA)))
      stop("Sensitivity and specifity must be fix values!")
    }
  } 

}

################################################################################
################################################################################

#-----------------------------------------------------------------------------
# input for rrisk.BayesZIP
#-----------------------------------------------------------------------------

checkInputZIP <- function(data,
                          prior.lambda, prior.pi,
                          chains, burn, thin, update,
                          workdir, plots
                          )
{
  
  checkMCparams(chains, burn, thin, update)
  checkEnvironParams(workdir, plots)
  
  if (missing(data) | missing(prior.lambda) | missing(prior.pi))
  { on.exit(return(invisible(NA)))
    stop("INVALID INPUT, missing one or more of the function arguments: 'data', 
         'prior.lambda', 'prior.pi'!", call. = FALSE)
  }
  
  if (!is.numeric(data) |!is.numeric(prior.lambda) | !is.numeric(prior.pi))
  { on.exit(return(invisible(NA)))
    stop("INVALID INPUT, one or more of the following arguments is not numeric: 
         'data', 'prior.lambda', 'prior.pi'!",
         call. = FALSE)
  }
  
  if (length(prior.lambda) != 2 | length(prior.pi) != 2)
  { on.exit(return(invisible(NA)))
    stop("INVALID INPUT, the arguments 'prior.lambda' and 'prior.pi' should be 
         of length 2!", call. = FALSE)
  }
 
  #-----------------------------------------------------------------------------
  # plausibility check on data
  #-----------------------------------------------------------------------------
  if (length(data) < 10)
  { on.exit(return(invisible(NA)))
    stop("Data set too small for this purpose!", call. = FALSE)
  }
  
  if (min(data) < 0)
  { on.exit(return(invisible(NA)))
    stop("Negative counts in data set are not allowed!", call. = FALSE)
  }
  
  if (any(abs(round(data) - data) != 0))
  { on.exit(return(invisible(NA)))
    stop("Data set contains non-integer values!", call. = FALSE)
  }
  
  #-----------------------------------------------------------------------------
  # plausibility check on priors
  #-----------------------------------------------------------------------------
  if (min(prior.pi) <= 0)
  { on.exit(return(invisible(NA)))
    stop("Parameters of the beta prior distribution for prevalence should not be
         negative!", call. = FALSE)
  }
  
  if (min(prior.lambda) < 0)
  { on.exit(return(invisible(NA)))
    stop("Parameters of the uniform prior distribution used as prior for the 
         Poisson parameter should be strictly positive!", call. = FALSE)
  }
}

################################################################################
################################################################################

#-----------------------------------------------------------------------------
# input for rrisk.BayesZINB
#-----------------------------------------------------------------------------
checkInputZINB <-  function(data,
                            prior.pi,
                            chains, burn, thin, update,
                            workdir, plots
                            )
{
  checkMCparams(chains, burn, thin, update)
  checkEnvironParams(workdir, plots)
  
  if (missing(data) | missing(prior.pi))
  { on.exit(return(invisible(NA)))
    stop("INVALID INPUT, missing one or more of the function arguments: 'data', 
         'prior.pi'!", call. = FALSE)
  }
  
  if (!is.numeric(data) | !is.numeric(prior.pi))
  { on.exit(return(invisible(NA)))
    stop("INVALID INPUT, one or more of the following arguments is not numeric: 
         'data', 'prior.pi'!",
         call. = FALSE)
  }
  
  if (length(prior.pi) != 2)
  { on.exit(return(invisible(NA)))
    stop("INVALID INPUT, 'prior.pi' should be of length 2!", call. = FALSE)
  }
  
  #-----------------------------------------------------------------------------
  # plausibility check on data
  #-----------------------------------------------------------------------------
  if (length(data) < 10)
  { on.exit(return(invisible(NA)))
    stop("Data set too small for this purpose!", call. = FALSE)
  }
  
  if (min(data) < 0)
  { on.exit(return(invisible(NA)))
    stop("Negative counts in data set are not allowed!", call. = FALSE)
  }
  
  if (any(abs(round(data) - data) != 0))
  { on.exit(return(invisible(NA)))
    stop("Data set contains non-integer values!", call. = FALSE)
  }
  
  #-----------------------------------------------------------------------------
  # plausibility check on priors
  #-----------------------------------------------------------------------------
  if (min(prior.pi) <= 0)
  { on.exit(return(invisible(NA)))
    stop("Parameters of the beta prior distribution for prevalence cannot be negative!", call. = FALSE)
  }
  
}
