#library(rriskBayes2)
test_rrisk.BayesPEM <- function(misclass, gui=FALSE){
  switch(misclass,
         "individual"={
           n <- 100
           x <- 14
           k <- 1
           pi_prior <- c(1, 1)
           se_prior <- c(64, 4)
           sp_prior <- c(94, 16)
           misclass <- "individual"
         },
         "individual-fix-sp"={
           n <- 100
           x <- 14
           k <- 1
           pi_prior <- c(1, 1)
           se_prior <- c(64, 4)
           sp_prior <- 1
           misclass <- "individual-fix-sp"
         },
         "individual-fix-se"={
           n <- 100
           x <- 14
           k <- 1
           pi_prior <- c(1, 1)
           sp_prior <- c(24, 5)
           se_prior <- 0.8
           misclass <- "individual-fix-se"
         },
         "individual-fix-se-sp"={
           n <- 100
           x <- 14
           k <- 1
           pi_prior <- c(1, 1)
           sp_prior <- 1
           se_prior <- 0.8
           misclass <- "individual-fix-se-sp"
         },
         "pool"={
           n <- 100
           k <- 10
           x <- 34
           pi_prior <- c(1, 1)
           se_prior <- c(64, 4)
           sp_prior <- c(94, 16)
           misclass <- "pool"
         },   
         "pool-fix-se"={
           n <- 100
           x <- 14
           k <- 1
           pi_prior <- c(1, 1)
           sp_prior <- c(24, 5)
           se_prior <- 0.8
           misclass <- "pool-fix-se"
         },
         "pool-fix-sp"={
           n <- 100
           x <- 14
           k <- 4
           pi_prior <- c(1, 1)
           sp_prior <- 1
           se_prior <- c(64, 4)
           misclass <- "pool-fix-sp"
         },
         "pool-fix-se-sp"={
           n <- 100
           x <- 14
           k <- 4
           pi_prior <- c(1, 1)
           sp_prior <- 1
           se_prior <- 0.7
           misclass <- "pool-fix-se-sp"
         },
         "compare"={
           n <- 100
           x <- 14
           k <- 4
           pi_prior <- c(1, 1)
           se_prior <- c(64, 4)
           sp_prior <- c(94, 16)
           misclass <- "compare"
         })

  if(!gui)
  # Estimate using Bayes model at individual level
   resPEM <- rrisk.BayesPEM(x = x, n = n, k = k, 
                            prior.se = se_prior,
                            prior.sp = sp_prior,
                            prior.pi = pi_prior,
                            misclass = misclass, 
                            chains=3,
                            burn = 4000,
                            thin = 1,
                            update = 10000,
                            workdir=getwd(),
                            plot=TRUE
                            )
  else
    resPEM <- PEMGUI(x = x, n = n, k = k, 
                     prior.se = se_prior,
                     prior.sp = sp_prior,
                     prior.pi = pi_prior,
                     chains=3,#
                     burn = 4000,
                     thin = 1,
                     update = 10000
                     )
  
  return(resPEM)
}


################################################################################
################################################################################
# generate ZIP data
# Beispieldaten generieren:


test_rrisk.BayesZIP <- function(gui=FALSE) {
  set.seed(42)
  
  n_true_neg <- 200
  n_true_pos <- 300
  n <- n_true_pos + n_true_neg
  lambda_true <- 0.5
  
  y_neg <- rep(0, n_true_neg)
  y_pos <- rpois(n_true_pos, lambda_true)
  
  # Data:
  y <- c(y_pos, y_neg)
  
  # Priors:
  lambda_prior <- c(0, 100)
  pi_prior     <- c(1, 1)
  
  lambda_prior =c(0, 100)
  pi_prior=c(1, 1)
  
  if(!gui)
  resZIP <- rrisk.BayesZIP(data = y,
                           prior.lambda = lambda_prior,
                           prior.pi = pi_prior,
                           simulation = FALSE,
                           chains = 3,
                           burn = 4000,
                           thin = 1,
                           update = 10000,
                           workdir = getwd(),
                           plots = TRUE
  )
  else
    resZIP <- ZIPGUI(data=y, prior.lambda = lambda_prior,
           prior.pi = pi_prior,
           chains = 3,
           burn = 4000,
           thin = 1,
           update = 10000)
  
  return(resZIP)
}


test_rrisk.BayesZINB <- function(gui = FALSE){
  # Beispieldaten generieren:
  #set.seed(42)
  
  n_true_neg <- 60
  n_true_pos <- 33
  n <- n_true_pos + n_true_neg
  
  prev_true <- n_true_pos / n
  
  a <- 6
  b <- 2
  lambda_true <- rgamma(n_true_pos, a, b)
  
  y_neg <- rep(0, n_true_neg)
  y_pos <- rpois(n_true_pos, lambda_true)
  y <- c(y_pos, y_neg)
  
  pi_prior     <- c(1, 1)
  
  
  if(!gui)
    resZINB <- rrisk.BayesZINB(data = y,
                               prior.pi = pi_prior,
                               simulation = FALSE,
                               chains = 3,
                               burn = 4000,
                               thin = 1,
                               update = 10000,
                               workdir = getwd(),
                               plots = TRUE
    )
  
  return(resZINB)
  
}






