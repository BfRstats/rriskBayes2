
# convergence check -------------------------------------------------------

checkPSRF <- function(x){
  ##check if psrf is smaller than target psrf
  res <- any(x$psrf$psrf[,"Point est."]< x$psrf$psrf.target)
  return(res)
}



