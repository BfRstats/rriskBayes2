
# convergence check -------------------------------------------------------

checkPSRF <- function(x){
  ##check if psrf is smaller than target psrf
  res <- any(x$psrf$psrf[,"Point est."]< x$psrf$psrf.target)
  return(res)
}

# binomial variable for zero-inflated models ------------------------------

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

