


################################################################################
################################################################################
#' @name plotDiag
#' @aliases plotDiag
#' @title Auxiliary function
#' @usage plotDiag(nodes, plotnumber)
#' @description Auxiliary function function for plotting graphical convergence
#'  diagnostics for Bayes models
#' @usage diagnostics(nodes, plots=FALSE)
#' @param S TO DESCRIBE
#' @param nodes character string, the name of parameters(s)
#' @param plotnumber single numerical value (1 for plotting Gelman-Rubin convergence
#'  statistic (\code{samplesBgr()}), 2 for plotting history (\code{plotDensity}),
#'  3 for plotting autocorrelation plots (\code{plotAutoC})) and 4 for plotting density plots
#' @keywords manip
#' @export

plotDiag <- function(S, plotnumber=1){
  switch(plotnumber,
         #Gelman-Rubin - Old: rrisk.sampleBgr
         "1" = ggs_Rhat(S) ,

         #Traceplot - Old: plotHistory
         "2" = ggs_traceplot(S),


         #Autocorrelation - Old: plotAutoC
         "3" = ggs_autocorrelation(S),

         #Density - Old: plotDensity
         "4" = ggs_density(S)
  )
}


