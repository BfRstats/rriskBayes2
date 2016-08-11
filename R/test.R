#library(rriskBayes2)
test_rrisk.BayesPEM <- function(misclass){
  n <- 100
  k <- 5
  x <- 56
  dat <- list(n = n, k = k, x = x)
  data_type <- ifelse(k == 1, "individual", "pooled")
  adapt <- 1000
  startburnin <- 4000
  burn <- adapt+startburnin
  
  #' # Priors
  #' Priors are used in the model string to define prior beta distributions.
  #' case individual
  switch(misclass,
         "individual"={
           pi_prior <- c(1, 1)
           se_prior <- c(86,16)
           sp_prior <- c(24, 5)
           misclass <- "individual"
         },
         "individual-fix-sp"={
           pi_prior <- c(1, 1)
           se_prior <- c(86,16)
           sp_prior <- 0.7
           misclass <- "individual-fix-sp"
         },
         "individual-fix-se"={
           pi_prior <- c(1, 1)
           sp_prior <- c(24, 5)
           se_prior <- 0.8
           misclass <- "individual-fix-se"
         },
         "pool"={
           pi_prior <- c(1, 1)
           se_prior <- c(64, 4)
           sp_prior <- c(94, 16)
           misclass <- "pool"
         },
         "compare"={
           pi_prior <- c(1, 1)
           se_prior <- c(64, 4)
           sp_prior <- c(94, 16)
           misclass <- "compare"
         })
 print(pi_prior)
 print(se_prior)
 print(sp_prior)
 print(misclass)
  # Estimate using Bayes model at individual level
  resPEM1 <- rrisk.BayesPEM(x=x, n=n, k=k, 
                            prior.se=se_prior,#
                            prior.sp=sp_prior,#
                            prior.pi=pi_prior,#
                            misclass=misclass, #
                            chains=3,#
                            burn=burn,
                            thin = 3,
                            update=10000,#
                            workdir=getwd(),#
                            plot=FALSE#
  )
}
################################################################################
################################################################################
# generate ZIP data
# Beispieldaten generieren:
set.seed(42)

n_true_neg <- 200
n_true_pos <- 300
n <- n_true_pos + n_true_neg
lambda_true <- 0.5

y_neg <- rep(0, n_true_neg)
y_pos <- rpois(n_true_pos, lambda_true)
y <- c(y_pos, y_neg)

# Daten:
jags_data <- list(y = y, n = n)
# PriorS:
lambda_prior <- c(0, 100)
pi_prior     <- c(1, 1)


resZIP <- rrisk.BayesZIP(data = y, prior.lambda = lambda_prior, prior.pi = pi_prior)

# estimate using Bayes model for zero inflated data
resZIP@results

prev_true <- n_true_pos / n
n_false0 <- sum(y_pos == 0)
n_pos <- sum(y_pos != 0)
# 
# zinb.data <- rep(0, n)
# zinb.data[sample(1:n, n*pi, replace=FALSE)] <- rnbinom(n*pi, size=65, prob=0.6)
# resZINB <- rrisk.BayesZINB(data = zinb.data)
# 



