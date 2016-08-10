
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
#' case individual
pi_prior <- c(1, 1)
se_prior <- c(86,16)
sp_prior <- c(24, 5)
misclass <- "individual"

#' case fixed sp
pi_prior <- c(1, 1)
se_prior <- c(86,16)
sp_prior <- 0.7
misclass <- "individual-fix-sp"

#' case fixed se
pi_prior <- c(1, 1)
sp_prior <- c(24, 5)
se_prior <- 0.8
misclass <- "individual-fix-se"

#' case pooled
pi_prior <- c(1, 1)
se_prior <- c(64, 4)
sp_prior <- c(94, 16)
misclass <- "pool"

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

################################################################################
################################################################################
# generate ZIP data
pi <- 0.8 #0.01
n <- 200
lambda <- 3.5
zip.data <- rep(0, n)
zip.data[sample(1:n, n*pi, replace=FALSE)] <- rpois(n*pi, lambda=lambda)
resZIP <- rrisk.BayesZIP(data = zip.data)
# estimate using Bayes model for zero inflated data
resZIP@results

rnbinom(n, size, prob, mu)

zinb.data <- rep(0, n)
zinb.data[sample(1:n, n*pi, replace=FALSE)] <- rnbinom(n*pi, size=65, prob=0.6)
resZINB <- rrisk.BayesZINB(data = zinb.data)

# 
# 
# n <- 10
# lambda  <-  runif(1, 0, 100)
# pi  <-   0.8 #dbeta(1, 1)
# library(MASS)
# parms <- fitdistr(zip.data, "poisson")
# y <- mu <- I <- rep(0, n)
# mu[1] <- parms$estimate
# 
# check <- function() {
#   for (i in 1:n) {
#     mu[i] <- I[i] * lambda
#     y[i]  <- rpois(n = 1, lambda = mu[i])
#     I[i] <- rbinom(1, 1, pi)
#   }
#   return(list(I=I,mu=mu,y=y))
# }


