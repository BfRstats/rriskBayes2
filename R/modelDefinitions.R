# defining the models for rrisk.BayesPEM ---------------------------------------

modelDefinitionPEM <- function(misclass) {
  
  #individual
  if (misclass == "individual") {
    model_string <- function(pi_prior, se_prior, sp_prior) {
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
    
    #individual fixed specifity
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
  # fixed sensitivity 
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
    
    #fixed sensitivity and fixed specifity
  } else if (misclass == "individual-fix-se-sp") {
    model_string <- function(pi_prior, se_fix, sp_fix) {
      sprintf(
        "model {
        pi ~  dbeta(%g, %g)                  # prevalence
        se <- %g                             # fixed sensitivity
        sp <-  %g                            # specificity
        ap <- pi*se + (1-pi)*(1-sp)          # apparent prevalence
        x  ~  dbin(ap, n)                    # number of positive samples
        
        #inits# pi
        #monitor# pi
        }",
        
        pi_prior[1],
        pi_prior[2],
        se_fix,
        sp_fix
      )
    }
    #pooled samples
  } else if (misclass == "pool") {
      model_string <- function(pi_prior, seP_prior, spP_prior) {
        sprintf(
        "model {
        pi  ~  dbeta(%g, %g)                 # prevalence
        se ~  dbeta(%g, %g)                  # (pooled test) sensitivity
        sp ~  dbeta(%g, %g)                  # (pooled test) specificity
        p.neg <- pow(1-pi, k)                # probability of a disease-free pool
        ap  <- (1-p.neg)*se + p.neg*(1-sp)   # apparent prevalence
        x   ~  dbin(ap, n)                   # number of positive pools
          
        #inits# pi, se, sp
        #monitor# pi, se, sp
        }",
          
        pi_prior[1],
        pi_prior[2],
        seP_prior[1],
        seP_prior[2],
        spP_prior[1],
        spP_prior[2]
      )
      }
      #pooled samples with fixed sensitivity
  } else if (misclass == "pool-fix-se") {
    model_string <- function(pi_prior, seP_fix, spP_prior) {
      sprintf(
        "model {
        pi  ~  dbeta(%g, %g)                 # prevalence
        se <-   %g                           # (pooled test) sensitivity
        sp ~  dbeta(%g, %g)                  # (pooled test) specificity
        p.neg <- pow(1-pi, k)                # probability of a disease-free pool
        ap  <- (1-p.neg)*se + p.neg*(1-sp)   # apparent prevalence
        x   ~  dbin(ap, n)                   # number of positive pools
          
        #inits# pi, sp
        #monitor# pi, sp
        }",
        
        pi_prior[1],
        pi_prior[2],
        seP_fix,
        spP_prior[1],
        spP_prior[2]
      )
    }
    #pooled samples with fixed specifity
  } else if (misclass == "pool-fix-sp") {
    model_string <- function(pi_prior, seP_prior, spP_fix) {
      sprintf(
        "model {
        pi  ~  dbeta(%g, %g)                 # prevalence
        se  ~  dbeta(%g, %g)                 # (pooled test) sensitivity
        sp <-  %g                            # (pooled test) specificity
        p.neg <- pow(1-pi, k)                # probability of a disease-free pool
        ap  <- (1-p.neg)*se + p.neg*(1-sp)   # apparent prevalence
        x   ~  dbin(ap, n)                   # number of positive pools
        
        #inits# pi, se
        #monitor# pi, se
    }",
        
        pi_prior[1],
        pi_prior[2],
        seP_prior[1],
        seP_prior[2],
        spP_fix
      )
    }
    #pooled samples with fixed sensitivity and fixed specifity
  } else if (misclass == "pool-fix-se-sp") {
    model_string <- function(pi_prior, seP_fix, spP_fix) {
      sprintf(
        "model {
        pi  ~  dbeta(%g, %g)                 # prevalence
        se  <- %g                            # (pooled test) sensitivity
        sp  <- %g                            # (pooled test) specificity
        p.neg <- pow(1-pi, k)                # probability of a disease-free pool
        ap  <- (1-p.neg)*se + p.neg*(1-sp)   # apparent prevalence
        x   ~  dbin(ap, n)                   # number of positive pools
        
        #inits# pi
        #monitor# pi
    }",
        
        pi_prior[1],
        pi_prior[2],
        seP_fix,
        spP_fix
      )
    }
    #compare 
  } else if (misclass == "compare"){
        model_string <- function(pi_prior, se_prior, sp_prior) {
        sprintf(
        "model {
        pi1 ~ dbeta(%g, %g)
        pi2 ~ dbeta(%g, %g)
        se ~ dbeta(%g, %g)
        sp ~ dbeta(%g, %g)
            
        p.neg <- pow(1-pi1, k)
        p.pos <- (1-p.neg)*se + p.neg*(1-sp)
        x1 ~ dbin(p.pos, n)
            
        ap <- pi2*se + (1-pi2)*(1-sp)
        p <- 1- pow(1-ap,k)
        x2 ~ dbin(p,n)
            
        d <- pi1 - pi2
            
        #inits# pi1, pi2, se, sp
        #monitor# pi1, pi2, se, sp, d
        }"
      ,
        pi_prior[1],
        pi_prior[2],
        pi_prior[1],
        pi_prior[2],
        se_prior[1],
        se_prior[2],
        sp_prior[1],
        sp_prior[2]
      )
   }
  }
  return(model_string)
}


# defining the model for rrisk.BayesZIP -----------------------------------

modelDefinitionZIP <- function(pi_prior, lambda_prior) {
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


# defining the model for rrisk.BayesZINB ----------------------------------

modelDefinitionZINB <- function(pi_prior) {
  sprintf(
        "model {
        pi ~ dbeta(%g, %g)
        dam ~ dgamma(0.01,0.01)
        db ~ dgamma(0.01,0.01)
    
        for (i in 1:n) {
        y[i]  ~ dpois(mu[i])
        mu[i] <- I[i] * lambda[i]
        I[i] ~ dbern(pi)
        lambda[i] ~ dgamma(dam, db)
        }
    
        #monitor# pi, dam, db
        }", 
        
        pi_prior[1], 
        pi_prior[2]
      )
   }
