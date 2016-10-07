
# plotDiag ----------------------------------------------------------------

#' @name plotDiag
#' @aliases plotDiag
#' @title Auxiliary function
#' @usage plotDiag(x, plotnumber)
#' @description Auxiliary function function for plotting graphical convergence
#'  diagnostics for Bayes models
#' @usage diagnostics(x, plots = FALSE)
#' @param x runjags object
#' @param plotnumber single numerical value (1 for plotting Gelman-Rubin convergence
#'  statistic (\code{samplesBgr()}), 2 for plotting history (\code{plotDensity}),
#'  3 for plotting autocorrelation plots (\code{plotAutoC})) and 4 for plotting density plots
#' @keywords manip
#' @export

plotDiag <- function(x, plotnumber = "all", ...){
  nvars <- length(x$monitor)
  
  switch(plotnumber,
         #Gelman-Rubin and Density - Old: plotnumber 1
         "1" = plot(x, plot.type = "trace", ...), #try(gelman.plot(x)),
         
         #Traceplot - Old: plotHistory plotnumber 2
         "2" = plot(x, plot.type = "ecdf", ...),

         #Autocorrelation - Old: plotAutoC  plotnumber 3
         "3" = plot(x, plot.type = "hist", ...),

         #Density
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
#' \code{\link{rrisk.BayesZIP}} and \code{\link{rrisk.BayesPEM}}.
#'
#' @details The argument \code{nodes} denotes the node(s) to be used for diagnostic plots 
#' of the MCMC chains. The user is interactively requested to confirm whether 
#' the convergence has been reached. In this case the function returns 
#' \code{TRUE} otherwise \code{FALSE}.
#' \cr \cr
#' This function is not intended to be called directly but is internally called
#' within \code{\link{rrisk.BayesZIP}} or \code{\link{rrisk.BayesPEM}}.
#'
#' @name diagnostics
#' @aliases diagnostics
#' @title Diagnostic plots for MCMC models provided by rrisk.BayesPEM and rrisk.BayesZIP functions
#' @usage diagnostics(nodes, plots = FALSE)
#' @param nodes character string, the name of parameters(s)
#' @param plots logical, if \code{TRUE} the diagnostic plots will be displayed in separate windows
#' @return Returns \code{TRUE} if the user confirms convergence. Otherwise the 
#' function returns \code{FALSE}.
# @note Some notes...
#' @seealso For more details see documentation to the functions \code{\link{samplesBgr}}, 
#' \code{\link{plotHistory}} and \code{\link{plotDensity}} from the package \pkg{BRugs}.
#' @keywords manip
#' @export
# @references Literatur to be added...
#' @examples
#' \dontrun{diagnostics(nodes)}

diagnostics<-function(x, plots = FALSE)
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
  allFrame<-tkframe(diagPlotWindow)
  headingFont2<-tkfont.create(weight = "bold",size = 10)
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
  
  
  # output
  
  return(get("output",envir = tempEnvir))
} #  end of function diagnostics()



