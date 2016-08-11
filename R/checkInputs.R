#-----------------------------------------------------------------------------
# checking input for rrisk.BayesPEM
#-----------------------------------------------------------------------------
checkInputPEM <- function(x,
                          n,
                          k,
                          prior.pi,
                          prior.se,
                          prior.sp,
                          chains,
                          burn,
                          thin,
                          update,
                          misclass,
                          workdir,
                          plots
) {
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
  
  if(chains<=0 | burn<=0 | update<=0 | thin <1)
  { on.exit(return(invisible(NA)))
    stop("INVALID INPUT, one or more of the following arguments is not positive: 'chains', 'burn', 'update', 'thin'!",call.=FALSE)
  }
  
  if(!is.element(misclass,c("pool","individual","individual-fix-sp", "individual-fix-se")))
  { on.exit(return(invisible(NA)))
    stop("INVALID INPUT, misclass should be 'pool','individual', 'individual-fix-sp' or individual-fix-se!",call.=FALSE)
  }
  
  if (missing(x) | missing(n))
  { on.exit(return(invisible(NA)))
    stop("INVALID INPUT, missing one or more of the following arguments: 'x', 'n'!", call. = FALSE)
  }
  
  if(!is.numeric(chains) | !is.numeric(burn) | !is.numeric(update) | !is.numeric(x) | !is.numeric(n))
  { on.exit(return(invisible(NA)))
    stop("INVALID INPUT, one or more of the following arguments is not numeric: 'chains', 'burn', 'update', 'x', 'n' !",call.=FALSE)
  }
  
  if(length(prior.sp)+length(prior.se)+length(prior.pi) < 5)
  { on.exit(return(invisible(NA)))
    stop("Only one variable ('se' or 'sp') can be fix!",call.=FALSE)
  }
  
  if (misclass == "pool") {
    if(missing(k) | !is.numeric(k))
    { on.exit(return(invisible(NA)))
      stop("INVALID INPUT, missing or non-numeric 'k'!", call. = FALSE)
    }
    
    max.l <- max(length(x),length(n),length(k)) # variable pool size disabled
    if((length(x) == max.l)*(length(n) == max.l)*(length(k) == max.l) != 1)
    { on.exit(return(invisible(NA)))
      stop("INVALID INPUT, equal length of arguments 'x', 'n' and 'k' required.",call.=FALSE)
    }
    
    if (missing(prior.pi) | missing(prior.se) | missing(prior.sp))
    { on.exit(return(invisible(NA)))
      stop("INVALID INPUT, missing one or more of the following prior arguments: 'prior.pi', 'prior.se', 'prior.sp'!", call. = FALSE)
    }
    
    if(length(prior.pi)!=2 | length(prior.se)!=2 | length(prior.sp)!=2)
    { on.exit(return(invisible(NA)))
      stop("Two parameters for the beta prior distribution for prevalence, sensitivity and specificity are required",call.=FALSE)
    }
    
  } else if (misclass == "individual-fix-sp") {
    if (length(prior.sp) != 1) 
    { on.exit(return(invisible(NA)))
      stop("Specifity must be a fix value!")
    }
    
    
  } else if (misclass == "individual-fix-se") {
    if (length(prior.se) != 1) 
    { on.exit(return(invisible(NA)))
      stop("Sensitivity must be a fix value!")
    }
    
  } else if (misclass == "individual" & !is.null(prior.pi))
  {
    warning("Your beta values are ignored. Beta prior is automatically set to c(1,1)!")
  }
  
  if(min(prior.pi)<=0 | min(prior.se)<=0 | min(prior.sp)<=0)
  { on.exit(return(invisible(NA)))
    stop("Parameters of the beta prior distribution for prevalence, sensitivity and specificity should be strictly positive!",call.=FALSE)
  }
}

################################################################################
################################################################################

#-----------------------------------------------------------------------------
# checking input for rrisk.BayesZIP
#-----------------------------------------------------------------------------
checkInputZIP <- function(data,
                          prior.lambda,
                          prior.pi,
                          simulation,
                          chains,
                          burn,
                          thin,
                          update,
                          workdir,
                          plots){
  if (missing(data))
  { on.exit(return(invisible(NA)))
    stop("INVALID INPUT, missing one or more of the function arguments: 'data', 'prior.lambda', 'prior.pi'!", call. = FALSE)
  }
  
  if (!is.numeric(data) |!is.numeric(prior.lambda) | !is.numeric(prior.pi) | !is.numeric(chains) | !is.numeric(burn) | !is.numeric(update))
  { on.exit(return(invisible(NA)))
    stop("INVALID INPUT, one or more of the following arguments is not numeric: 'data', 'prior.lambda', 'prior.pi', 'chains', 'burn', 'update'!",
         call. = FALSE)
  }
  
  if (!is.logical(plots))
  { on.exit(return(invisible(NA)))
    stop("INVALID INPUT, the argument 'plots' should be of type logical!", call. = FALSE)
  }
  
  if (chains <= 0 | burn <= 0 | update <= 0 | thin < 1)
  { on.exit(return(invisible(NA)))
    stop("INVALID INPUT, one or more of the following arguments is not positive: 'chains', 'burn', 'update', 'thin'!", call. = FALSE)
  }
  
  if (length(prior.lambda) != 2 | length(prior.pi) != 2)
  { on.exit(return(invisible(NA)))
    stop("INVALID INPUT, the arguments 'prior.lambda' and 'prior.pi' should be of length 2!", call. = FALSE)
  }
  
  try.setwd <- try(setwd(workdir), silent = TRUE)
  if (inherits(try.setwd, "try-error"))
  { on.exit(return(invisible(NA)))
    error.mess <- paste("INVALID INPUT, the working directory could not be found! Your input is \n", workdir)
    stop(error.mess, call. = FALSE)
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
    stop("Parameters of the beta prior distribution for prevalence should not be negative!", call. = FALSE)
  }
  
  if (min(prior.lambda) < 0)
  { on.exit(return(invisible(NA)))
    stop("Parameters of the uniform prior distribution used as prior for the Poisson parameter should be strictly positive!", call. = FALSE)
  }
}
