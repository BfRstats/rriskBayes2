checkInput <- function(x,
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
    {
      on.exit(return(invisible(NA)))
      stop("INVALID INPUT, missing or non-numeric 'k'!", call. = FALSE)
    }

    max.l <- max(length(x),length(n),length(k)) # variable pool size disabled
    if((length(x) == max.l)*(length(n) == max.l)*(length(k) == max.l) != 1)
    { on.exit(return(invisible(NA)))
      stop("INVALID INPUT, equal length of arguments 'x', 'n' and 'k' required.",call.=FALSE)
    }

    if (missing(prior.pi) | missing(prior.se) | missing(prior.sp))
    {
      on.exit(return(invisible(NA)))
      stop("INVALID INPUT, missing one or more of the following prior arguments: 'prior.pi', 'prior.se', 'prior.sp'!", call. = FALSE)
    }

    if(length(prior.pi)!=2 | length(prior.se)!=2 | length(prior.sp)!=2)
    { on.exit(return(invisible(NA)))
      stop("Two parameters for the beta prior distribution for prevalence, sensitivity and specificity are required",call.=FALSE)
    }


  } else if (misclass == "individual-fix-sp") {
    if (length(prior.sp) != 1) {
      on.exit(return(invisible(NA)))
      stop("Specifity must be a fix value!")
    }


  } else if (misclass == "individual-fix-se") {
    if (length(prior.se) != 1) {
      on.exit(return(invisible(NA)))
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

inits_function <-
  function(chain, misclass) {
    # max number of chains: 5
    
    pi <- c(0.2, 0.5, 0.8, 0.9, 0.7)[chain]
    se <- c(0.9, 0.8, 0.7, 0.6, 0.75)[chain]
    sp <- c(0.7, 0.9, 0.8, 0.85, 0.95)[chain]
    
    .RNG.seed <- c(1, 2, 3, 4, 5)[chain]
    .RNG.name <- c(
      "lecuyer::RngStream",
      "base::Super-Duper",
      "base::Wichmann-Hill",
      "base::Mersenne-Twister",
      "base::Marsaglia-Multicarry"
    )[chain]
    
    if (misclass == "individual-fix-se") {
      inits.list <- list(
        .RNG.seed = .RNG.seed,
        .RNG.name = .RNG.name,
        pi = pi,
        sp = sp
      )
    } else if (misclass == "individual-fix-sp") {
      inits.list <- list(
        .RNG.seed = .RNG.seed,
        .RNG.name = .RNG.name,
        pi = pi,
        se = se
      )
      
    } else if (misclass == "pool" | misclass == "individual") {
      inits.list <- list(
        .RNG.seed = .RNG.seed,
        .RNG.name = .RNG.name,
        pi = pi,
        seP = se,
        spP = sp
      )
    }
    
    return(inits.list)
  }


################################################################################
################################################################################

#-----------------------------------------------------------------------------
# write model
#-----------------------------------------------------------------------------
writeModel <- function(misclass) {
  if (misclass == "individual") {
   model <- 
      "model{
      pi ~ dbeta(1, 1)
      
      se ~ dbeta(prior.se[1],prior.se[2])
      
      sp ~ dbeta(prior.sp[1],prior.sp[2])
      
      ap <- pi*se + (1-pi)*(1-sp)
      
      x ~ dbin(ap,n)
  }"
    
} else if (misclass == "individual-fix-sp") {
  model <- 
    "model{
    pi ~ dbeta(prior.pi[1],prior.pi[2])
    
    se ~ dbeta(prior.sp[1],prior.sp[2])
    
    sp ~ fix
    
    ap <- pi*se + (1-pi)*(1-sp)
    
    x ~ dbin(p,n)
}"
    
  } else if (misclass == "individual-fix-se") {
    model <- 
      "model{
      pi ~ dbeta(prior.pi[1],prior.pi[2])
      
      se ~ fix
      
      sp ~ dbeta(prior.sp[1],prior.sp[2])
      
      ap <- pi*se + (1-pi)*(1-sp)
      
      x ~ dbin(p,n)
}"
    
    } else if (misclass == "compare") {
      model <-
        "model{
        pi1 ~ dbeta(prior.pi[1],prior.pi[2])
        
        pi2 ~ dbeta(prior.pi[1],prior.pi[2])
        
        se ~ dbeta(prior.se[1],prior.se[2])
        
        sp ~ dbeta(prior.sp[1],prior.sp[2])
        
        x1 <- x
        
        x2 <- x
        
        p.neg <- pow(1-pi1,k)
        
        p.pos <- (1-p.neg)*se + p.neg*(1-sp)
        
        x1 ~ dbin(p.pos,n)
        
        ap <- pi2*se + (1-pi2)*(1-sp)
        
        p <- 1- pow(1-ap,k)
        
        x2 ~ dbin(p,n)
}"
      
      } else if (misclass == "pool") {
        model <- 
          "model{
          
          pi ~ dbeta(prior.pi[1],prior.pi[2])
          
          se ~ dbeta(prior.se[1],prior.se[2])
          
          sp ~ dbeta(prior.sp[1],prior.sp[2])
          
          p.neg <- pow(1-pi,k)
          
          p <- (1-p.neg)*se + p.neg*(1-sp)
          
          x ~ dbin(p,n)
}"
        
      }

  return(model)
  
  }

