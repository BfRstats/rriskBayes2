
# common for all models ---------------------------------------------------
# select random number generators (only up to 5 chains) 
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

# individual for the models ---------------------------------------------------
# binomial variable for zero-inflated models (ZIP and  ZINB)-------------------

inits_chainZI <- function(chain, data, seed_nr = 33) {
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

# defining the initialization parameters for rrisk.BayesPEM ---------------

inits_functionPEM <- function(chain, misclass) {
  # max number of chains: 5
  
  pi <- c(0.2, 0.5, 0.8, 0.9, 0.7)[chain]
  se <- c(0.9, 0.8, 0.7, 0.6, 0.75)[chain]
  sp <- c(0.7, 0.9, 0.8, 0.85, 0.95)[chain]
  
  randNr <- RNGs(chain)
  .RNG.seed <- randNr$.RNG.seed
  .RNG.name <- randNr$.RNG.name
  
  if (misclass == "individual" | misclass == "pool") {
    inits.list <- list(
      .RNG.seed = .RNG.seed,
      .RNG.name = .RNG.name,
      pi = pi,
      se = se,
      sp = sp
    )
  } else if(misclass == "individual-fix-se" | misclass == "pool-fix-se") {
    inits.list <- list(
      .RNG.seed = .RNG.seed,
      .RNG.name = .RNG.name,
      pi = pi,
      sp = sp
    )
  } else if (misclass == "individual-fix-sp"  | misclass == "pool-fix-sp") {
    inits.list <- list(
      .RNG.seed = .RNG.seed,
      .RNG.name = .RNG.name,
      pi = pi,
      se = se
    )
  } else if (misclass == "individual-fix-se-sp" | misclass == "pool-fix-se-sp") {
    inits.list <- list(
      .RNG.seed = .RNG.seed,
      .RNG.name = .RNG.name,
      pi = pi
    )
  } else if (misclass == "compare") {
    inits.list <- list(
      .RNG.seed = .RNG.seed,
      .RNG.name = .RNG.name,
      pi1 = pi,
      pi2 = pi,
      se = se,
      sp = sp
    )
  } 
  return(inits.list)
}


# defining the initialization parameters for rrisk.BayesZIP ---------------

inits_functionZIP <-  function(chain, data) {
  # max number of chains: 5
  pi <- c(0.2, 0.5, 0.8, 0.9, 0.7)[chain]
  lambda <- c(1, 80, 0.2, 2, 70)[chain]
  I <- inits_chainZI(chain = chain, data = data)
  
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


# defining the initialization parameters for rrisk.BayesZINB --------------

inits_functionZINB <-  function(chain, data) {
  # max number of chains: 5
  pi <- c(0.2, 0.5, 0.8, 0.9, 0.7)[chain]
  dam <- c(1, 80, 0.2, 2, 70)[chain]
  db <- c(2, 3, 10, 20, 2)[chain]
  
  I <- inits_chainZI(chain = chain, data = data)
  
  randNr <- RNGs(chain)
  .RNG.seed <- randNr$.RNG.seed
  .RNG.name <- randNr$.RNG.name
  
  inits.list <- list(
    .RNG.seed = .RNG.seed,
    .RNG.name = .RNG.name,
    pi = pi,
    dam = dam,
    db = db, 
    I = I
  )
  return(inits.list)
}
