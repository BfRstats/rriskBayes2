# write model rrisk.BayesPEM ----------------------------------------------


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
    
    sp <-  fix
    
    ap <- pi*se + (1-pi)*(1-sp)
    
    x ~ dbin(ap,n)
    }"
  } else if (misclass == "individual-fix-se") {
    model <-
      "model{
    pi ~ dbeta(prior.pi[1],prior.pi[2])
    
    se <-  fix
    
    sp ~ dbeta(prior.sp[1],prior.sp[2])
    
    ap <- pi*se + (1-pi)*(1-sp)
    
    x ~ dbin(ap,n)
    }"
  } else if (misclass == "individual-fix-se-sp") {
    model <-
      "model{
    pi ~ dbeta(prior.pi[1],prior.pi[2])
    
    se <-  fix
    
    sp <-  fix
    
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
      
    p.neg <- pow(1-pi1,k)
      
    p.pos <- (1-p.neg)*se + p.neg*(1-sp)
      
    x1 ~ dbin(p.pos,n)
      
    ap <- pi2*se + (1-pi2)*(1-sp)
      
    p <- 1- pow(1-ap,k)
      
    x2 ~ dbin(p,n)
      
    d <- pi1 - pi2
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
  } else if (misclass == "pool-fix-se") {
      model <-
      "model{
    pi ~ dbeta(prior.pi[1],prior.pi[2])
    
    se <- fix
        
    sp ~ dbeta(prior.sp[1],prior.sp[2])
        
    p.neg <- pow(1-pi,k)
        
    p <- (1-p.neg)*se + p.neg*(1-sp)
        
    x ~ dbin(p,n)
    }"
  } else if (misclass == "pool-fix-sp") {
      model <-
      "model{
    pi ~ dbeta(prior.pi[1], prior.pi[2])
        
    se <- dbeta(prior.se[1],prior.se[2])
        
    sp <- fix
        
    p.neg <- pow(1-pi,k)
        
    p <- (1-p.neg)*se + p.neg*(1-sp)
        
    x ~ dbin(p,n)
    }"
  } else if (misclass == "pool-fix-se-sp") {
      model <-
      "model{
    pi ~ dbeta(prior.pi[1],prior.pi[2])
    
    se <- fix
    
    sp <- fix
    
    p.neg <- pow(1-pi,k)
    
    p <- (1-p.neg)*se + p.neg*(1-sp)
    
    x ~ dbin(p,n)
  }"
 }
  return(model)
}

# write model rrisk.BayesZIP ----------------------------------------------


writeModelZIP <- function() {
  model <-
    "model{
    pi ~ dbeta(prior.pi[1], prior.pi[2])

    lambda ~ dunif(prior.lambda[1], prior.lambda[2])
  
    for (i in 1:n) {
      y[i]  ~ dpois(mu[i])
      mu[i] <- I[i] * lambda
      I[i] ~ dbern(pi)
    }
  }"
  return(model)
  }


# write model rrisk.BayesZINB ---------------------------------------------


writeModelZINB <- function() {
  model <-
    "model{
    pi  ~ dbeta(prior.pi[1], prior.pi[2])

    dam ~ dgamma(0.01,0.01)

    db  ~ dgamma(0.01,0.01)
  
    for (i in 1:n) {
      y[i]  ~ dpois(mu[i])
      mu[i] <- I[i] * lambda[i]
      I[i] ~ dbern(pi)
      lambda[i] ~ dgamma(dam,db)
    }
  }"
  return(model)
}
