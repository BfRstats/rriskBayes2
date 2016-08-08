#library(rriskBayes2)
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
pi_prior <- c(1, 1)
se_prior <- c(86,16)
sp_prior <- c(24, 5)

# Estimate using Bayes model at individual level
resPEM1 <- rrisk.BayesPEM(x=x, n=n,
                          prior.se=se_prior,#
                          prior.sp=sp_prior,#
                          prior.pi=pi_prior,#
                          misclass="individual", #
                          chains=3,#
                          burn=burn,#
                          update=10000,#
                          workdir=getwd(),#
                          plot=FALSE#
)
