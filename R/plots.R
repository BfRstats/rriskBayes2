################################################################################
################################################################################
#' @name plotDiag
#' @aliases plotDiag
#' @title Auxiliary function
#' @usage plotDiag(nodes, plotnumber)
#' @description Auxiliary function function for plotting graphical convergence
#'  diagnostics for Bayes models
#' @usage diagnostics(nodes, plots=FALSE)
#' @param x runjags object
#' @param nodes character string, the name of parameters(s)
#' @param plotnumber single numerical value (1 for plotting Gelman-Rubin convergence
#'  statistic (\code{samplesBgr()}), 2 for plotting history (\code{plotDensity}),
#'  3 for plotting autocorrelation plots (\code{plotAutoC})) and 4 for plotting density plots
#' @keywords manip
#' @export

plotDiag <- function(x, plotnumber = "all", ...){
  nvars <- length(x$monitor)
  
  switch(plotnumber,
         #Gelman-Rubin and Density - Old: plotnumber 1
         "1" = try(gelman.plot(x)),
         
         #Traceplot - Old: plotHistory plotnumber 2
         "2" = plot(x, plot.type = "trace", ...),

         #Autocorrelation - Old: plotAutoC  plotnumber 3
         "3" = plot(x, plot.type = "autocorr", ...),

         #Density
         "4" = plot(x, plot.type = "density", ...),

         "all" = {
           try(gelman.plot(x))
           plot(x, plot.type="density", ...)
           plot(x, plot.type= "trace", ...)
           plot(x, plot.type = "autocorr", ...)
         }

  )
}


