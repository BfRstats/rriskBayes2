#library(rriskBayes2)
test_rrisk.BayesPEM <- function(misclass){
  switch(misclass,
         "individual"={
           n <- 100
           x <- 14
           k <- NULL
           pi_prior <- c(1, 1)
           se_prior <- c(64, 4)
           sp_prior <- c(94, 16)
           misclass <- "individual"
         },
         "individual-fix-sp"={
           n <- 100
           x <- 14
           k <- NULL
           pi_prior <- c(1, 1)
           se_prior <- c(64, 4)
           sp_prior <- 1
           misclass <- "individual-fix-sp"
         },
         "individual-fix-se"={
           n <- 100
           x <- 14
           k <- NULL
           pi_prior <- c(1, 1)
           sp_prior <- c(24, 5)
           se_prior <- 0.8
           misclass <- "individual-fix-se"
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
         "compare"={
           n <- 100
           x <- 14
           k <- NULL
           pi_prior <- c(1, 1)
           se_prior <- c(64, 4)
           sp_prior <- c(94, 16)
           misclass <- "compare"
         })

  dat <- list(n = n, k = k, x = x)
  data_type <- ifelse(k == 1, "individual", "pooled")

  # Estimate using Bayes model at individual level
   resPEM <- rrisk.BayesPEM(x = x, n = n, k = k, 
                            prior.se = se_prior,#
                            prior.sp = sp_prior,#
                            prior.pi = pi_prior,#
                            misclass = misclass, #
                            chains=3,#
                            burn = 4000,
                            thin = 1,
                            update = 10000,#
                            workdir=getwd(),#
                            plot=TRUE#
  )
  return(resPEM)
}
################################################################################
################################################################################
# generate ZIP data
# Beispieldaten generieren:


test_rrisk.BayesZIP <- function() {
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
  
  resZIP <- rrisk.BayesZIP( data = y,
                            prior.lambda = lambda_prior,
                            prior.pi = pi_prior,
                            simulation = FALSE,
                            chains = 3,
                            burn = 4000,
                            thin = 1,
                            update = 10000,
                            workdir = getwd(),
                            plots = FALSE
 )
  
  return(resZIP)
}





