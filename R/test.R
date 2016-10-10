#library(rriskBayes2)
# test risk.BayesPEM ------------------------------------------------------


#' @title Quick test of rrisk.BayesPEM
#' @description The function \code{test_rrisk.BayesPEM} is a quick test of the \code{\link{rrisk.BayesPEM}} function for all cases of \code{misclass} without setting parameters. x, n and k are sampled, all other parameters are fixed. 
#' @param misclass See \code{\link{misclass}} character with legal character entries: \cr
#' \code{individual}, \code{individual-fix-se}, \code{individual-fix-sp}, \code{individual-fix-se-sp}, \cr
#' \code{pool}, \code{pool-fix-se}, \code{pool-fix-sp} or \code{pool-fix-se-sp}.\cr
#' \code{fix-se}: fixed sensitivity\cr
#' \code{fix-sp}: fixed specifity\cr
#' \code{fix-se-sp}: fixed sensitivity AND fixed specifity
#' @return Returns an instance of the 
#' \code{\linkS4class{bayesmodelClass}} class. 
#' @export
#'
#' @examples
#' \donttest{
#'------------------------------------------
#' test examples
#'------------------------------------------
#' 
#' res <- test_rrisk.BayesPEM(misclass = "individual")
#' res <- test_rrisk.BayesPEM(misclass = "compare")
#' }
test_rrisk.BayesPEM <- function(misclass){
  switch(misclass,
         "individual"={
           x <- sample(0:100, 1)
           n <- sample(x:400, 1)
           k <- 1
           pi_prior <- c(1, 1)
           se_prior <- c(64, 4)
           sp_prior <- c(94, 16)
           misclass <- "individual"
         },
         "individual-fix-sp"={
           x <- sample(0:100, 1)
           n <- sample(x:400, 1)
           k <- 1
           pi_prior <- c(1, 1)
           se_prior <- c(64, 4)
           sp_prior <- 1
           misclass <- "individual-fix-sp"
         },
         "individual-fix-se"={
           x <- sample(0:100, 1)
           n <- sample(x:400, 1)
           k <- 1
           pi_prior <- c(1, 1)
           sp_prior <- c(24, 5)
           se_prior <- 0.8
           misclass <- "individual-fix-se"
         },
         "individual-fix-se-sp"={
           x <- sample(0:100, 1)
           n <- sample(x:400, 1)
           k <- 1
           pi_prior <- c(1, 1)
           sp_prior <- 1
           se_prior <- 0.8
           misclass <- "individual-fix-se-sp"
         },
         "pool"={
           x <- sample(0:100, 1)
           n <- sample(x:400, 1)
           k <- sample(5:50, 1)
           pi_prior <- c(1, 1)
           se_prior <- c(64, 4)
           sp_prior <- c(94, 16)
           misclass <- "pool"
         },   
         "pool-fix-se"={
           x <- sample(0:100, 1)
           n <- sample(x:400, 1)
           k <- sample(5:50, 1)
           pi_prior <- c(1, 1)
           sp_prior <- c(24, 5)
           se_prior <- 0.8
           misclass <- "pool-fix-se"
         },
         "pool-fix-sp"={
           x <- sample(0:100, 1)
           n <- sample(x:400, 1)
           k <- sample(5:50, 1)
           pi_prior <- c(1, 1)
           sp_prior <- 1
           se_prior <- c(64, 4)
           misclass <- "pool-fix-sp"
         },
         "pool-fix-se-sp"={
           x <- sample(0:100, 1)
           n <- sample(x:400, 1)
           k <- sample(5:50, 1)
           pi_prior <- c(1, 1)
           sp_prior <- 1
           se_prior <- 0.7
           misclass <- "pool-fix-se-sp"
         },
         "compare"={
           x <- sample(0:100, 1)
           n <- sample(x:400, 1)
           k <- sample(5:50, 1)
           pi_prior <- c(1, 1)
           se_prior <- c(64, 4)
           sp_prior <- c(94, 16)
           misclass <- "compare"
         })

  # Estimate using Bayes model at individual level
   resPEM <- rrisk.BayesPEM(x = x, n = n, k = k, 
                            prior.se = se_prior,
                            prior.sp = sp_prior,
                            prior.pi = pi_prior,
                            misclass = misclass, 
                            chains=3,
                            burn = 4000,
                            thin = 3,
                            plots=FALSE
                            )
  
  return(resPEM)
}


# test rrisk.BayesZIP -----------------------------------------------------

#' @title Quick test of \code{rrisk.BayesZIP} 
#' @description The function \code{\link{test_rrisk.BayesZIP}} is a quick test of the rrisk.BayesZIP function without setting parameters. Number of true negatives and true positives are sampled, all other parameters were fixed (\code{prior.lambda = c(0, 100)}, \code{prior.pi = c(1,1)}) . 
#' @return Returns an instance of the 
#' \code{\linkS4class{bayesmodelClass}} class. 
#' @export
#' @examples 
#' \donttest{
#'------------------------------------------
#' test example
#'------------------------------------------
#' 
#' res <- test_rrisk.BayesZIP()
#' }
test_rrisk.BayesZIP <- function() {
  
  # generate ZIP data
  n_true_neg <- sample(1:500, 1) 
  n_true_pos <- sample(1:500, 1) 
  n <- n_true_pos + n_true_neg
  lambda_true <- 0.5
  
  y_neg <- rep(0, n_true_neg)
  y_pos <- rpois(n_true_pos, lambda_true)
  
  y <- c(y_pos, y_neg)
  
  # Priors:
  lambda_prior <- c(0, 100)
  pi_prior     <- c(1, 1)

  resZIP <- rrisk.BayesZIP(data = y,
                           prior.lambda = lambda_prior,
                           prior.pi = pi_prior,
                           simulation = FALSE,
                           chains = 3,
                           burn = 4000,
                           thin = 1,
                           update = 10000,
                           plots = FALSE)

  return(resZIP)
}


# test rrisk.BayesZINB ----------------------------------------------------

#' @title Quick test of the \code{rrisk.BayesZINB} function
#' @description The function is a quick test of the \code{\link{rrisk.BayesZINB}} function without setting parameters. Number of true negatives (n_true_neg) and true positives (n_true_pos) are sampled, all other parameters are fixed \cr
#' \code{prior.pi = c(1,1)} \cr
#' \code{lambda_true <- rgamma(n_true_pos, a=6, b=2)} \cr
#' \code{y_neg <- rep(0, n_true_neg)} \cr
#' \code{y_pos <- rpois(n_true_pos, lambda_true)} \cr
#' \code{y <- c(y_pos, y_neg)}
#' @return Returns an instance of the 
#' \code{\linkS4class{bayesmodelClass}} class. 
#' @export
#' @examples 
#' \donttest{
#'------------------------------------------
#' test example
#'------------------------------------------
#' 
#' res <- test_rrisk.BayesZINB()
#' }
test_rrisk.BayesZINB <- function(){

  n_true_neg <- sample(1:500, 1) 
  n_true_pos <- sample(1:500, 1) 
  n <- n_true_pos + n_true_neg
  
  a <- 6; b <- 2
  lambda_true <- rgamma(n_true_pos, a, b)
  pi_prior     <- c(1, 1)
  
  y_neg <- rep(0, n_true_neg)
  y_pos <- rpois(n_true_pos, lambda_true)
  y <- c(y_pos, y_neg)
  
  
  return(resZINB)
}






