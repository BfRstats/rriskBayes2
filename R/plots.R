


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

plotDiag <- function(x, plot.type="all"){
  switch(plot.type,
         #Gelman-Rubin - Old: rrisk.sampleBgr
         "gelman" = gelman.plot(x) ,

         #Traceplot - Old: plotHistory
         "trace" = plot(x, plot.type= "trace"),

         #Autocorrelation - Old: plotAutoC
         "autocorr" = plot(x, plot.type = "autocorr"),

         #Histogram
         "histogram" = plot(x, plot.type = "histogram"),

         #Density - Old: plotDensity
         "density" = {
           mcmcobj <- suppressWarnings(as.mcmc(x))
           par(mfrow=c(1,3))
           densplot(mcmcobj)
           par(mfrow=c(1,1))
         },

         "all" =
         {
           par(mfrow = c(2,3))
           gelman.plot(x, auto.layout = FALSE)
           mcmcobj <- suppressWarnings(as.mcmc(x))
           densplot(mcmcobj)
           par(mfrow=c(1,1))
           plot(x, plot.type= "trace")
           plot(x, plot.type = "autocorr")
          # plot(x, plot.type = "histogram")
         }

  )
}


