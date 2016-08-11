################################################################################
################################################################################

#-----------------------------------------------------------------------------
# functions applicable for all methods
#-----------------------------------------------------------------------------

################################################################################
################################################################################
#-----------------------------------------------------------------------------
# select random number generators (only up to 5 chains)
#-----------------------------------------------------------------------------
RNGs <- function(chain)
{
  .RNG.seed <- c(1, 2, 3, 4, 5)[chain]
  .RNG.name <- c(
    "lecuyer::RngStream",
    "base::Super-Duper",
    "base::Wichmann-Hill",
    "base::Mersenne-Twister",
    "base::Marsaglia-Multicarry"
  )[chain]
  
  return(list(.RNG.seed = .RNG.seed, .RNG.name = .RNG.name))
}

################################################################################
################################################################################

#-----------------------------------------------------------------------------
# functions applicable for specific methods
#-----------------------------------------------------------------------------

################################################################################
################################################################################

#-----------------------------------------------------------------------------
# defining the models for rrisk.BayesPEM
#-----------------------------------------------------------------------------

modelFunctionPEM <- function(misclass) {
  if (misclass == "individual") {
    model_string <- function(pi_prior = pi_prior, se_prior, sp_prior) {
      sprintf(
        "model {
        pi ~  dbeta(%g, %g)                 # prevalence
        se ~  dbeta(%g, %g)                 # sensitivity
        sp ~  dbeta(%g, %g)                 # specificity
        ap <- pi*se + (1-pi)*(1-sp)         # apparent prevalence
        x  ~  dbin(ap, n)                   # number of positive samples
        
        #inits# pi, se, sp
        #monitor# pi, se, sp
      }",
        pi_prior[1],
        pi_prior[2],
        se_prior[1],
        se_prior[2],
        sp_prior[1],
        sp_prior[2]
      )
    }
} else if (misclass == "individual-fix-sp") {
  model_string <- function(pi_prior, se_prior, sp_fix) {
    sprintf(
      "model {
        pi ~  dbeta(%g, %g)                  # prevalence
        se ~  dbeta(%g, %g)                  # sensitivity
        sp <- %g                             # fixed specificity
        ap <- pi*se + (1-pi)*(1-sp)          # apparent prevalence
        x  ~  dbin(ap, n)                    # number of positive samples
      
        #inits# pi, se
        #monitor# pi, se
      }",
        pi_prior[1],
        pi_prior[2],
        se_prior[1],
        se_prior[2],
        sp_fix
      )
    }
} else if (misclass == "individual-fix-se") {
    model_string <- function(pi_prior, se_fix, sp_prior) {
      sprintf(
        "model {
        pi ~  dbeta(%g, %g)                  # prevalence
        se <- %g                             # fixed sensitivity
        sp ~  dbeta(%g, %g)                  # specificity
        ap <- pi*se + (1-pi)*(1-sp)          # apparent prevalence
        x  ~  dbin(ap, n)                    # number of positive samples
        
        #inits# pi, sp
        #monitor# pi, sp
      }",
        pi_prior[1],
        pi_prior[2],
        se_fix,
        sp_prior[1],
        sp_prior[2]
      )
    }
} else if (misclass == "pool") {
      model_string <- function(pi_prior, seP_prior, spP_prior) {
        sprintf(
          "model {
        pi  ~  dbeta(%g, %g)                 # prevalence
        seP ~  dbeta(%g, %g)                 # (pooled test) sensitivity
        spP ~  dbeta(%g, %g)                 # (pooled test) specificity
        p.neg <- pow(1-pi, k)                # probability of a disease-free pool
        ap  <- (1-p.neg)*seP + p.neg*(1-spP) # probability for a pool
        # to test positive
        x   ~  dbin(ap, n)                   # number of positive pools
          
        #inits# pi, seP, spP
        #monitor# pi, seP, spP
      }",
        pi_prior[1],
        pi_prior[2],
        seP_prior[1],
        seP_prior[2],
        spP_prior[1],
        spP_prior[2]
      )
    }
  }
 return(model_string)
}

################################################################################
################################################################################

#-----------------------------------------------------------------------------
# defining the models for rrisk.BayesZIP
#-----------------------------------------------------------------------------

modelFunctionZIP <- function(pi_prior, lambda_prior) {
  sprintf(
    "model{
        pi  ~  dbeta(%g, %g)
        lambda  ~  dunif(%g, %g)
       
        for (i in 1:n){
          y[i]  ~ dpois(mu[i])
          mu[i] <- I[i] * lambda
          I[i] ~ dbern(pi)
        } 

        #monitor# pi, lambda
      }", 
        pi_prior[1],
        pi_prior[2],
        lambda_prior[1],
        lambda_prior[2]
        )
  }

################################################################################
################################################################################


#-----------------------------------------------------------------------------
# defining the models for rrisk.BayesZIP
#-----------------------------------------------------------------------------

modelFunctionZINB <- function(r_prior, p_prior, pi_prior) {
  sprintf(
    "model{
    r  ~  dunif(%g, %g)
    p  ~  dunif(%g, %g)
    pi  ~  dbeta(%g, %g)
    
    for (i in 1:n){
     y[i]  ~ dnegbin(mu[i], r)  
     mu[i] <- I[i] * p
     I[i] ~ dbern(pi)
    } 
    
    #inits# r, p, pi
    #monitor# r, p, pi
}", 
      r_prior[1],
      r_prior[2],
      p_prior[1],
      p_prior[2],
      pi_prior[1],
      pi_prior[2]
  )
  }


################################################################################
################################################################################

#-----------------------------------------------------------------------------
# defining the initialization parameters for rrisk.BayesPEM
#-----------------------------------------------------------------------------

inits_functionPEM <- function(chain, misclass) {
  # max number of chains: 5
  
  pi <- c(0.2, 0.5, 0.8, 0.9, 0.7)[chain]
  se <- c(0.9, 0.8, 0.7, 0.6, 0.75)[chain]
  sp <- c(0.7, 0.9, 0.8, 0.85, 0.95)[chain]
  
  randNr <- RNGs(chain)
  .RNG.seed <- randNr$.RNG.seed
  .RNG.name <- randNr$.RNG.name
  
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
  } else if (misclass == "pool" |
             misclass == "individual") {
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
# defining the initialization parameters for rrisk.BayesZIP
#-----------------------------------------------------------------------------
inits_chainZIP <- function(chain, data, seed_nr = 33) {
  set.seed(seed_nr)
  # inits für chain 1: nur Einsen
  I1 <- rep(1, data$n)
  # inits für chain 2: Einsen, wenn y > 0 ist, sonst Nullen
  I2 <- as.numeric(data$y != 0)
  # inits für chains 3-5: zufälliger 0/1-Vektor, aber alle Positionen,
  # für die y > 0 ist, werden auf Eins gesetzt
  idx_pos <- which(data$y != 0)
  I3 <- sample(c(0, 1), data$n, replace = TRUE)
  I3[idx_pos] <- 1
  I4 <- sample(c(0, 1), data$n, replace = TRUE)
  I4[idx_pos] <- 1
  I5 <- sample(c(0, 1), data$n, replace = TRUE)
  I5[idx_pos] <- 1
  I <- list(I1, I2, I3, I4, I5)[[chain]]
  return(I)
}

inits_functionZIP <-  function(chain, data) {
  # max number of chains: 5
  
  pi <- c(0.2, 0.5, 0.8, 0.9, 0.7)[chain]
  lambda <- c(1, 80, 0.2, 2, 70)[chain]
  I <- inits_chainZIP(chain = chain, data = data)
  
  randNr <- RNGs(chain)
  .RNG.seed <- randNr$.RNG.seed
  .RNG.name <- randNr$.RNG.name
  
  inits.list <- list(
    .RNG.seed = .RNG.seed,
    .RNG.name = .RNG.name,
    pi = pi,
    lambda = lambda,
    I = I
  )
  return(inits.list)
}

################################################################################
################################################################################

#-----------------------------------------------------------------------------
# defining the initialization parameters for rrisk.BayesZIP
#-----------------------------------------------------------------------------

inits_functionZINB <-  function(chain) {
  # max number of chains: 5
  
  r <- c(4, 5, 6, 7, 8)[chain]
  p <- c(0.9, 0.8, 0.7, 0.6, 0.75)[chain]
  pi <- c(0.2, 0.5, 0.8, 0.9, 0.7)[chain]
 
  randNr <- RNGs(chain)
  .RNG.seed <- randNr$.RNG.seed
  .RNG.name <- randNr$.RNG.name
  
  inits.list <- list(
    .RNG.seed = .RNG.seed,
    .RNG.name = .RNG.name,
    r = r,
    p = p,
    pi = pi
  )
  return(inits.list)
}

################################################################################
################################################################################

#-----------------------------------------------------------------------------
# write model rrisk.BayesPEM
#-----------------------------------------------------------------------------

writeModelPEM <- function(misclass) {
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
  
        x ~ dbin(ap,n)
      }"
  } else if (misclass == "individual-fix-se") {
    model <-
      "model{
        pi ~ dbeta(prior.pi[1],prior.pi[2])
    
        se ~ fix
    
        sp ~ dbeta(prior.sp[1],prior.sp[2])
    
        ap <- pi*se + (1-pi)*(1-sp)
    
        x ~ dbin(ap,n)
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

################################################################################
################################################################################

#-----------------------------------------------------------------------------
# write model rrisk.BayesZIP
#-----------------------------------------------------------------------------

writeModelZIP <- function() {
  model <-
      "model{
        pi ~ dbeta(prior.pe[1], prior.pi[2])
        lambda ~ dunif(prior.lambda[1], prior.lambda[2])

        for (i in 1:n) {
        y[i]  ~ dpois(mu[i])
        mu[i] <- I[i] * lambda
        I[i] ~ dbern(pi)
        }
        }"
  return(model)
}

