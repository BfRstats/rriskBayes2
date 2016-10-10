
# plotDiag ----------------------------------------------------------------

#' @name plotDiag
#' @aliases plotDiag
#' @title Auxiliary function
#' @usage plotDiag(x, plotnumber)
#' @description Auxiliary function function for plotting graphical convergence
#'  diagnostics for Bayes models
#' @usage plotDiag <- function(x, plotnumber = "all", ...)
#' @param x Runjags object
#' @param plotnumber Single numerical value \cr
#' 1 for plotting the trace \cr
#' 2 for plotting the cumulative distribution function \cr
#' 3 for plotting the histogram \cr
#' 4 for plotting the autocorrelation
#' 5 for plotting all plots and additionally the Gelman and Rubin plot\cr
#' @param ... Further plotting parameters from the \code{plot} function of \code{\link{runjags}}. 
#' @keywords manip
#' @examples
#' \dontrun{
#' resPEM <- rrisk.BayesPEM(x = 14, n = 100, k = 4, prior.se = c(64, 4), prior.sp = c(94, 16))
#' plotDiag(resPEM@results, plotnumber = 3)
#' plotDiag(resPEM@results, plotnumber = "all")
#' }
#' @export

plotDiag <- function(x, plotnumber = "all", ...){
  nvars <- length(x$monitor)
  
  switch(plotnumber,
         #Traceplot - Old: plotHistory plotnumber 1
         "1" = plot(x, plot.type = "trace", ...), 
         
         #ECDF - Old: plotDensity  
         "2" = plot(x, plot.type = "ecdf", ...),

         #Histogram - Old: plotAutoC  plotnumber 2
         "3" = plot(x, plot.type = "hist", ...),

         #Autocorrelation - Old: plotAutoC  plotnumber 3
         "4" = plot(x, plot.type = "autocorr", ...),

         "all" = {
           try(gelman.plot(x))
           plot(x, plot.type="hist", ...)
           plot(x, plot.type="ecdf", ...)
           plot(x, plot.type= "trace", ...)
           plot(x, plot.type = "autocorr", ...)
         }
  )
}


# diagnostics -------------------------------------------------------------


#' @description This function provides a GUI for diagnostic plots to check convergence in Markov
#' chain Monte-Carlo (MCMC) models provided by 
#' \code{\link{rrisk.BayesZIP}}, \code{\link{rrisk.BayesZINB}} and \code{\link{rrisk.BayesPEM}}.
#' @details The argument \code{x} denotes the result of the MCMC models. The user is interactively requested to confirm whether 
#' the convergence has been reached. In this case the function returns 
#' \code{TRUE} otherwise \code{FALSE}.
#' \cr
#' The function is not intended to be called directly but is internally called
#' within the models.
#' @name diagnostics
#' @aliases diagnostics
#' @title Diagnostic plots for MCMC models provided by rrisk.BayesPEM, rrisk.BayesZIP and rrisk.BayesZINB functions
#' @usage diagnostics(x, plots = FALSE)
#' @param x Result slot returned by \code{\link{rrisk.BayesPEM}}, \code{\link{rrisk.BayesZIP}} or \code{\link{rrisk.BayesZINB}} of class \code{\link{runjags}}.
#' @param plots Logical, if \code{TRUE} the diagnostic plots will be displayed in separate windows
#' @export
#' @return Returns \code{TRUE} if the user confirms convergence. Otherwise the 
#' function returns \code{FALSE}.
#' @keywords manip
# @references Literatur to be added...
#' @examples
#' \dontrun{
#' resPEM <- rrisk.BayesPEM(x = 14, n = 100, k = 4, prior.se = c(64, 4), prior.sp = c(94, 16))
#' diagnostics(resPEM@results)
#' }

diagnostics <- function(x, plots = FALSE)
{ 
  # what happends by pressing "not converged" button
  onNO<-function(...)
  {
    assign("output", value = FALSE, envir = tempEnvir)
    tkdestroy(diagPlotWindow)
  } # end of onNO()
 
  
  # what happends by pressing "converged" button
  
  onYES <- function(...)
  {
    assign("output", value = TRUE, envir = tempEnvir)
    tkdestroy(diagPlotWindow)
  } # end of onNO()
  
  # display diagnostic plots in separate windows
  
  if(plots)
  {
    X11();plotDiag(x, plotnumber = 1)
    X11();plotDiag(x, plotnumber = 2)
    X11();plotDiag(x, plotnumber = 3)
  }
  
  # define help variables for output
  
  assign("tempEnvir",value = new.env(), envir = .GlobalEnv)
  assign("output", value = FALSE, envir = tempEnvir)
  assign("plotNumber", value = 1, envir = tempEnvir)
  
  diagPlotWindow<-tktoplevel(width = 100,height = 100)
  tkwm.title(diagPlotWindow,"Graphical convergence diagnostics for Bayesian models")
  tkwm.resizable(diagPlotWindow,0,0)  # fixed size, not resizeable
  tkwm.maxsize(diagPlotWindow,800,600)
  tkwm.minsize(diagPlotWindow,800,600)
  allFrame <- tkframe(diagPlotWindow)
  headingFont2<-tkfont.create(weight = "bold", size = 10)
  imgFrame<-tkframe(allFrame)

  #"trace", "ecdf", "histogram", "autocorr" plots
  vars <- x$monitor
  nvars <- length(vars)
  k <- 1
  for(i in 1:nvars) {
    assign(paste0("imgPlot", k), value = tkrplot(imgFrame,fun = function() 
      plot(x, vars = vars[k]), hscale = 1.6, vscale = 1.5), envir = tempEnvir)
    k <- k+1
  }
  #"crosscorrelation" plot
  if(nvars > 1)
    assign(paste0("imgPlot", nvars+1), value = tkrplot(imgFrame,fun = function() 
      plot(x, vars = vars, plot.type = "crosscorr"), hscale = 1.6, vscale = 1.5), 
      envir = tempEnvir)
  
  tkpack(imgFrame,fill = "both",expand = "yes")
  imgPlot1 <- get("imgPlot1", envir = tempEnvir)
  tkpack(imgPlot1)
  buttonsFrame<-tkframe(allFrame)
  nextButton<-ttkbutton(buttonsFrame,width = nchar("Show next plot")+2, 
                        text = "Show next plot", 
                        command = function(...) onNextplot(envir = tempEnvir, nvars = nvars))
  convButton <- ttkbutton(buttonsFrame,width = nchar("Converged")+2, 
                        text = "Converged",comman = onYES)
  notconvButton <- ttkbutton(buttonsFrame,width = nchar("Not converged")+2, 
                           text = "Not converged",command = onNO)
  tkpack(nextButton,convButton,notconvButton,padx = c(20,20),side = "left")
  tkpack(buttonsFrame,pady = c(0,0))
  tkpack(allFrame,fill = "both",expand = "yes",padx = c(15,15),pady = c(15,15))
  
  tkfocus(diagPlotWindow)
  tkwait.window(diagPlotWindow)
  
  return(get("output",envir = tempEnvir))
} #  end of function diagnostics()



