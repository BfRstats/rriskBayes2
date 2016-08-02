
################################################################################
################################################################################
#' @description This function calculates and plots the Gelman-Rubin convergence statistic,
#' as modified by Brooks and Gelman (1998).
#'
#' @details This function is an alias of the function \code{\link{samplesBgr}} from the
#' package \pkg{BRugs}. The original function was modified to create a diagnostic
#' plots in other format as it was implemented in the original version. For more
#' details see the function \code{\link{samplesBgr}} from the package \pkg{BRugs}.
#' \cr \cr
#' This function is not intended to be called directly but is internally called
#' by \code{\link{diagnostics}} within the functions
#' \code{\link{rrisk.BayesZIP}} and \code{\link{rrisk.BayesPEM}}.
#'
#' @name rrisk.samplesBgr
#' @aliases rrisk.samplesBgr
#' @title Plot the Gelman-Rubin convergence statistic
#' @usage rrisk.samplesBgr(node, beg=samplesGetBeg(), end=samplesGetEnd(),
#'  firstChain=samplesGetFirstChain(), lastChain=samplesGetLastChain(),
#'  thin=samplesGetThin(), bins=50, plot=TRUE, ...)
#' @param node character vector of length 1, name of a variable in the model
#' @param beg argument to select a slice of monitored values corresponding to iterations \code{beg:end}
#' @param end argument to select a slice of monitored values corresponding to iterations \code{beg:end}
#' @param firstChain argument to select a sub group of chains to calculate the Gelman-Rubin convergence statistics for. Number of chains must be larger than one
#' @param lastChain argument to select a sub group of chains to calculate the Gelman-Rubin convergence statistics for. Number of chains must be larger than one
#' @param thin only use every thin-th value of the stored sample for statistics
#' @param bins number of blocks
#' @param plot logical, whether to plot the BGR statistics or only return the values. If \code{TRUE}, values are returned invisibly
#' @param ... further graphical parameters as in \code{par}
# @return Returns ...
# @note Some notes...
# @seealso Nothing...
#' @keywords manip
#' @export
#' @examples
#' \dontrun{rrisk.samplesBgr("se")}

rrisk.samplesBgr<-function(node, beg=samplesGetBeg(), end=samplesGetEnd(),
                           firstChain=samplesGetFirstChain(), lastChain=samplesGetLastChain(),
                           thin=samplesGetThin(), bins=50, plot=TRUE, ...)
{
  #if (plot && is.null(ask))
  if (plot)
  {
    if (is.R())
      ask <- !((dev.cur() > 1) && !dev.interactive())
    else ask <- !((dev.cur() > 1) && !interactive())
  }
  oldBeg <- samplesGetBeg()
  oldEnd <- samplesGetEnd()
  oldFirstChain <- samplesGetFirstChain()
  oldLastChain <- samplesGetLastChain()
  oldThin <- samplesGetThin()
  on.exit({
    samplesSetBeg(oldBeg)
    samplesSetEnd(oldEnd)
    samplesSetFirstChain(oldFirstChain)
    samplesSetLastChain(oldLastChain)
    samplesSetThin(oldThin)})
  beg <- max(beg, modelAdaptivePhase())
  samplesSetBeg(beg)
  samplesSetEnd(end)
  samplesSetFirstChain(firstChain)
  samplesSetLastChain(lastChain)
  thin <- max(c(thin, 1))
  samplesSetThin(thin)
  mons <- samplesMonitors(node)
  result <- lapply(mons, plotBgr, bins = bins, plot = plot,...)
  names(result) <- mons
  if (plot)
    invisible(result)
  else return(result)
} # end of function rrisk.samplesBgr()


################################################################################
################################################################################
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
#' @usage diagnostics(nodes, plots=FALSE)
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

diagnostics<-function(nodes, plots=FALSE)
{
  #-----------------------------------------------------------------------------
  # what happends by pressing "not converged" button
  #-----------------------------------------------------------------------------
  onNO<-function(...)
  {
    assign("output",value=FALSE,envir=tempEnvir)
    tkdestroy(diagPlotWindow)
  } # end of onNO()
  #-----------------------------------------------------------------------------
  # what happends by pressing "converged" button
  #-----------------------------------------------------------------------------
  onYES<-function(...)
  {
    assign("output",value=TRUE,envir=tempEnvir)
    tkdestroy(diagPlotWindow)
  } # end of onNO()
  #-----------------------------------------------------------------------------
  # what happends by pressing "next plot" button
  #-----------------------------------------------------------------------------
  onNext<-function(...)
  { plotNumber<-get("plotNumber",envir=tempEnvir)
  if(plotNumber==1){
    tkpack.forget(imgPlot1)
    assign("plotNumber",value=2,envir=tempEnvir)
    tkpack(imgPlot2,side="top")
  } else if (plotNumber==2){
    tkpack.forget(imgPlot2)
    assign("plotNumber",value=3,envir=tempEnvir)
    tkpack(imgPlot3,side="top")
  } else if (plotNumber==3){
    tkpack.forget(imgPlot3)
    assign("plotNumber",value=1,envir=tempEnvir)
    tkpack(imgPlot1,side="top")
  }
  tkraise(diagPlotWindow)
  }  # end of onNext()

  #-----------------------------------------------------------------------------
  # display diagnostic plots in separate windows
  #-----------------------------------------------------------------------------
  if(plots)
  {
    X11();plotDiag(nodes,plotnumber=1)
    X11();plotDiag(nodes,plotnumber=2)
    X11();plotDiag(nodes,plotnumber=3)
  }

  #-----------------------------------------------------------------------------
  # define help variables for output
  #-----------------------------------------------------------------------------
  assign("tempEnvir",value=new.env(),envir=.GlobalEnv)
  assign("output",value=FALSE,envir=tempEnvir)
  assign("plotNumber",value=1,envir=tempEnvir)

  diagPlotWindow<-tktoplevel(width=100,height=100)
  tkwm.title(diagPlotWindow,"Graphical convergence diagnostics for Bayesian models")
  tkwm.resizable(diagPlotWindow,0,0)  # fixed size, not resizeable
  tkwm.maxsize(diagPlotWindow,800,600)
  tkwm.minsize(diagPlotWindow,800,600)
  allFrame<-tkframe(diagPlotWindow)
  headingFont2<-tkfont.create(weight="bold",size=10)
  imageFrame<-tkframe(allFrame)

  imgPlot1<-tkrplot(imageFrame,fun=function()plotDiag(nodes,plotnumber=1),hscale=2,vscale=1.4);
  imgPlot2<-tkrplot(imageFrame,fun=function()plotDiag(nodes,plotnumber=2),hscale=2.0,vscale=1.4)
  imgPlot3<-tkrplot(imageFrame,fun=function()plotDiag(nodes,plotnumber=3),hscale=1.0,vscale=1.4)

  tkpack(imageFrame,fill="both",expand="yes")
  tkpack(imgPlot1)
  buttonsFrame<-tkframe(allFrame)
  nextButton<-ttkbutton(buttonsFrame,width=nchar("Show next plot")+2, text="Show next plot",command=onNext)
  convButton<-ttkbutton(buttonsFrame,width=nchar("Converged")+2, text="Converged",comman=onYES)
  notconvButton<-ttkbutton(buttonsFrame,width=nchar("Not converged")+2, text="Not converged",command=onNO)
  tkpack(nextButton,convButton,notconvButton,padx=c(20,20),side="left")
  tkpack(buttonsFrame,pady=c(0,0))
  tkpack(allFrame,fill="both",expand="yes",padx=c(15,15),pady=c(15,15))

  tkfocus(diagPlotWindow)
  tkwait.window(diagPlotWindow)
  #-----------------------------------------------------------------------------
  # output
  #-----------------------------------------------------------------------------
  return(get("output",envir=tempEnvir))
} #  end of function diagnostics()


################################################################################
################################################################################
#' @name plotDiag
#' @aliases plotDiag
#' @title Auxiliary function
#' @usage plotDiag(nodes, plotnumber)
#' @description Auxiliary function function for plotting graphical convergence
#'  diagnostics for Bayes models
#' @usage diagnostics(nodes, plots=FALSE)
#' @param nodes character string, the name of parameters(s)
#' @param plotnumber single numerical value (1 for plotting Gelman-Rubin convergence
#'  statistic (\code{samplesBgr()}), 2 for plotting history (\code{plotDensity}) and
#'  3 for plotting autocorrelation plots (\code{plotAutoC}))
#' @keywords manip
#' @export

plotDiag<-function(nodes,plotnumber=1)
{ len<-length(nodes)
if(plotnumber==1){
  par(mfrow=c(2,len))
  for(i in 1:len){
    if(len<=3){
      maintext<-paste("Gelman-Rubin conv. \nstatistic for",nodes[i])
    } else {
      maintext<-paste("Gelman-Rubin conv. \nstatistic for",nodes[i])
    }
    suppressWarnings(rrisk.samplesBgr(nodes[i],main=maintext,cex.main=1))
    suppressWarnings(abline(h=1,lty="dashed"))
  }
  for(i in 1:len){
    suppressWarnings(plotDensity(nodes[i],main=paste("Density plot of ",nodes[i]),cex.main=1))
  }
} else if (plotnumber==2){
  if(len<=3){
    par(mfrow=c(len,1))
  } else  par(mfrow=c(2,3))
  #par(mfrow=c(len,1))
  for(i in 1:len){
    suppressWarnings(plotHistory(nodes[i],main=paste("Convergence monitored for node",nodes[i]),cex.main=1))
  }
} else if(plotnumber==3){
  if(len<=3){
    par(mfrow=c(len,1))
  } else  par(mfrow=c(2,3))
  for(i in 1:len){
    suppressWarnings(plotAutoC(nodes[i],main=paste("Autocorrelation plot for node",nodes[i]),cex.main=1))
  }
} # end if statement
} # end of function plotDiag()

################################################################################
################################################################################
#' @description This function provides a GUI for the function \link{rrisk.BayesZIP}.
#'
#' @name ZIPGUI
#' @aliases ZIPGUI
#' @title Bayes estimation of a zero inflated Poisson (ZIP) model
#' @usage ZIPGUI(data=NULL, prior.lambda=c(1,10), prior.pi=c(0.8,1),
#'  chains=3, burn=1000, update=10000, thin=1)
#' @param data a vector of numeric data, containing zeros, and of minimal length 10.
#' @param prior.lambda numeric vector containing minimum and maximum of a uniform
#' distribution used as prior for the Poisson parameter \code{lambda}, e.g. \cr \code{lambda} \eqn{\sim} \code{prior.lambda(*,*)=unif(*,*)}
#' @param prior.pi numeric vector containing parameters of a beta distribution
#' describing prior knowledge about prevalence (proportion of contaminated samples), e.g. \cr \code{pi} \eqn{\sim} \code{prior.pi(*,*)=beta(*,*)}
#' @param chains positive single numeric value, number of independent MCMC chains (default 3)
#' @param burn positive single numeric value, length of the burn-in period (default 1000)
#' @param update positive single numeric value, length of update iterations for estimation (default 10000)
#' @param thin positive single numeric value (default 1). The samples from every kth iteration will be used for
#'        inference, where k is the value of thin. Setting \code{thin > 1} can help to reduce the autocorrelation
#'        in the sample.
#' @return The function \code{ZIPGUI} returns an instance of the \code{\linkS4class{bayesmodelClass}}
#' class containing following information
#' \item{\code{convergence}}{logical, whether the model has converged (assessed by the user)}
#' \item{\code{results}}{data frame containing statitsics of the posterior distribution}
#' \item{\code{jointpost}}{data frame giving the joint posterior probability distribution}
#' \item{\code{nodes}}{names of the parameters jointly estimated by the Bayes model}
#' \item{\code{model}}{model in BRugs/Winbugs syntax as a character string}
#' \item{\code{chains}}{number of independent MCMC chains}
#' \item{\code{burn}}{length of burn-in period}
#' \item{\code{update}}{length of update iterations for estimation}
#' @note The convergence of the model will be entered by the user after the simulation process.
#' @seealso \code{\link{rrisk.BayesZIP}}
#' @keywords manip
#' @export
#' @references Bohning, D., E. Dietz, P. Schlattman, L. Mendonca, and U. Kirchner (1999).
#' The zero-inflated Poisson model and the decayed, missing and filled teeth index in
#' dental epidemiology. Journal of the Royal Statistical Society, Series A 162, 195-209.
#' @examples
#' \donttest{
#' data <- rpois(30, 4)
#' prior.lambda <- c(1, 10)
#' prior.pi <- c(0.8, 1)
#' ZIPGUI(data, prior.lambda, prior.pi)}


ZIPGUI <- function(data=NULL, prior.lambda=c(1, 10), prior.pi=c(0.8, 1),
                   chains=3, burn=1000, update=10000,thin=1){
  #-----------------------------------------------------------------------------
  # define buttons and their actions
  #-----------------------------------------------------------------------------
  # GUIDiag: function for plotting the two diagnosis plots
  GUIDiag<-function(nodes,plotnumber=1){
    len<-length(nodes)
    if(plotnumber==1){
      par(mfrow=c(2,len))
      for(i in 1:len){
        if(len<=3){
          maintext<-paste("Gelman-Rubin conv. \nstatistic for",nodes[i])
        } else {
          maintext<-paste("G.-R. statistic for",nodes[i])
        } # end if
        suppressWarnings(rrisk.samplesBgr(nodes[i],main=maintext,cex.main=1))
        suppressWarnings(abline(h=1,lty="dashed")  )
      } # end first for-loop
      for(i in 1:len){
        suppressWarnings(plotDensity(nodes[i],main=paste("Density plot of ",nodes[i]),cex.main=1))
      } # end second for-loop
    } else if (plotnumber==2){
      if(len<=3){
        par(mfrow=c(len,1))
      } else {
        par(mfrow=c(2,3))
      } #enf if-else
      for(i in 1:len){
        suppressWarnings(plotHistory(nodes[i],main=paste("Convergence monitored for node",nodes[i]),cex.main=1))
      } # emd for-loop
    } else if(plotnumber==3){
      if(len<=3){
        par(mfrow=c(len,1))
      } else {
        par(mfrow=c(2,3))
      } #enf if-else
      for(i in 1:len){
        suppressWarnings(plotAutoC(nodes[i],main=paste("Autocorrelation plot for node",nodes[i]),cex.main=0.5))
      } # emd for-loop
    }# end overall if
  } # end function GUIDiag()

  #-----------------------------------------------------------------------------
  # what happens by pressing RUN button
  #-----------------------------------------------------------------------------
  onRun <- function(...){
    prior.pi.l. <- as.numeric(tclvalue(tkget(prior.pi.l)))
    prior.pi.r. <- as.numeric(tclvalue(tkget(prior.pi.r)))
    prior.lambda.l. <- as.numeric(tclvalue(tkget(prior.lambda.l)))
    prior.lambda.r. <- as.numeric(tclvalue(tkget(prior.lambda.r)))
    chains. <- as.numeric(tclvalue(tkget(chainsEntry)))
    burn. <- as.numeric(tclvalue(tkget(burnEntry)))
    update. <- as.numeric(tclvalue(tkget(updateEntry)))
    thin. <- as.numeric(tclvalue(tkget(thinEntry)))

    ##try.result <- try(
    mod <- rrisk.BayesZIP(data=data, simulation=TRUE, prior.pi=c(prior.pi.l., prior.pi.r.),
                          prior.lambda=c(prior.lambda.l., prior.lambda.r.),
                          chains=chains., burn=burn., update=update.,thin=thin.)

    if(!is(mod,"bayesmodelClass")){
      assign("mod", value=mod, envir=envirZIP)
      tkdestroy(bayesZIPWindow)
      stop("Any error occured during fitting bayesian ZIP model",call.=FALSE)
    } else {
      ##,silent=TRUE)
      ##if (inherits(try.result, "try-error")) { cat("Error: The simulation could not be completed!\n")}
      assign("mod", value=mod, envir=envirZIP)

      if(exists("imgPlot1", where=envirZIP)){
        imgPlot1 <- get("imgPlot1", envir=envirZIP)
        tkpack.forget(imgPlot1)
      }
      if(exists("imgPlot2", where=envirZIP)){
        imgPlot2 <- get("imgPlot2", envir=envirZIP)
        tkpack.forget(imgPlot2)
      }
      if(exists("imgPlot3", where=envirZIP)){
        imgPlot2 <- get("imgPlot3", envir=envirZIP)
        tkpack.forget(imgPlot3)
      }

      nodes <- c("pi", "lambda")
      imgPlot1 <- tkrplot(imgFrame,fun=function()GUIDiag(nodes,plotnumber=1),hscale=1.4,vscale=1.4)
      imgPlot2 <- tkrplot(imgFrame,fun=function()GUIDiag(nodes,plotnumber=2),hscale=1.4,vscale=1.4)
      imgPlot3 <- tkrplot(imgFrame,fun=function()GUIDiag(nodes,plotnumber=3),hscale=1.4,vscale=1.4)

      assign("imgPlot1", value=imgPlot1, envir=envirZIP)
      assign("imgPlot2", value=imgPlot2, envir=envirZIP)
      assign("imgPlot3", value=imgPlot3, envir=envirZIP)

      tkpack(imgPlot1, side="top")
      assign("plotNumber",value=1,envir=envirZIP)

      tkraise(bayesZIPWindow)
    }
  } # end onRun() fucntion

  #-----------------------------------------------------------------------------
  # what happens by pressing RESET button
  #-----------------------------------------------------------------------------
  onReset <- function(...){
    tkconfigure(prior.pi.l, text=tclVar(prior.pi[1]))
    tkconfigure(prior.pi.r, text=tclVar(prior.pi[2]))
    tkconfigure(prior.lambda.l, text=tclVar(prior.lambda[1]))
    tkconfigure(prior.lambda.r, text=tclVar(prior.lambda[2]))
    tkconfigure(chainsEntry, text=tclVar(chains))
    tkconfigure(burnEntry, text=tclVar(burn))
    tkconfigure(updateEntry, text=tclVar(update))
    tkconfigure(thinEntry, text=tclVar(thin))
  } # end onReset() function

  #-----------------------------------------------------------------------------
  # what happens by pressing CANCEL button
  #-----------------------------------------------------------------------------
  onCancel <- function(...){
    tkdestroy(bayesZIPWindow)
  } # end onCancel() function

  #-----------------------------------------------------------------------------
  # what happens by pressing Nexplot button
  #-----------------------------------------------------------------------------
  onNextplot <- function(...){
    # get the current plot number (in 1,2) and switch to the other plot
    if(!exists("plotNumber", envirZIP)) return(FALSE) # if there was no model run yet

    plotNumber <- get("plotNumber", envir=envirZIP)
    if(plotNumber==1){
      imgPlot1 <- get("imgPlot1", envir=envirZIP)
      tkpack.forget(imgPlot1)
      assign("plotNumber",value=2,envir=envirZIP)
      imgPlot2 <- get("imgPlot2", envir=envirZIP)
      tkpack(imgPlot2, side="top")
    }  else if(plotNumber==2){
      imgPlot2 <- get("imgPlot2", envir=envirZIP)
      tkpack.forget(imgPlot2)
      assign("plotNumber",value=3,envir=envirZIP)
      imgPlot3 <- get("imgPlot3", envir=envirZIP)
      tkpack(imgPlot3, side="top")
    } else if(plotNumber==3){
      imgPlot3 <- get("imgPlot3", envir=envirZIP)
      tkpack.forget(imgPlot3)
      assign("plotNumber",value=1,envir=envirZIP)
      imgPlot1 <- get("imgPlot1", envir=envirZIP)
      tkpack(imgPlot1, side="top")
    }
    tkraise(bayesZIPWindow)
  } # end onNextplot() function

  #-----------------------------------------------------------------------------
  # what happens by pressing CONVERGENCE button
  #-----------------------------------------------------------------------------
  onConv <- function(...){
    mod <- get("mod", envir=envirZIP)
    mod@convergence <- TRUE
    assign("mod", value=mod, envir=envirZIP)
    tkdestroy(bayesZIPWindow)
  } # end fucntion onConv()

  #-----------------------------------------------------------------------------
  # what happens by pressing CONVERGENCE button
  #-----------------------------------------------------------------------------
  onNotconv <- function(...){
    mod <- get("mod", envir=envirZIP)
    mod@convergence <- FALSE
    assign("mod", value=mod, envir=envirZIP)
    tkdestroy(bayesZIPWindow)
  } # end fucntion onNotconv()

  #-----------------------------------------------------------------------------
  # define help varriable(s)
  #-----------------------------------------------------------------------------
  assign("envirZIP", value=new.env(), envir=.GlobalEnv)

  #-----------------------------------------------------------------------------
  # define GUI window and frames
  #-----------------------------------------------------------------------------
  bayesZIPWindow <- tktoplevel()
  tkwm.title(bayesZIPWindow, "GUI for the function rrisk.bayesZIP")
  tkwm.resizable(bayesZIPWindow, 0, 0) # fixed size
  tkwm.maxsize(bayesZIPWindow,800,500)
  tkwm.minsize(bayesZIPWindow,800,500 )

  mod <- new("bayesmodelClass")

  leftFrame <- tkframe(bayesZIPWindow)
  rightFrame <- tkframe(bayesZIPWindow)
  imgFrame <- tkframe(rightFrame,height=500,width=200)
  inputFrame <- tkframe(leftFrame)
  lButtonFrame <- tkframe(leftFrame)
  rButtonFrame <- tkframe(leftFrame)
  #rButtonFrame <- tkframe(rightFrame)

  #-----------------------------------------------------------------------------
  # define input fields
  #-----------------------------------------------------------------------------
  prior.lambda.l <- tkentry(inputFrame, text=tclVar(prior.lambda[1]), width=6)
  prior.lambda.lLabel <- tklabel(inputFrame, text="prior.lambda (unif)")
  prior.lambda.r <- tkentry(inputFrame, text=tclVar(prior.lambda[2]), width=6)

  prior.pi.l <- tkentry(inputFrame, text=tclVar(prior.pi[1]), width=6)
  prior.pi.lLabel <- tklabel(inputFrame, text="prior.pi (beta)")
  prior.pi.r <- tkentry(inputFrame, text=tclVar(prior.pi[2]), width=6)

  chainsEntry <- tkentry(inputFrame, text=tclVar(chains))
  chainsLabel <- tklabel(inputFrame, text="chains")

  burnEntry <- tkentry(inputFrame, text=tclVar(burn))
  burnLabel <- tklabel(inputFrame, text="burn")

  updateEntry <- tkentry(inputFrame, text=tclVar(update))
  updateLabel <- tklabel(inputFrame, text="update")

  thinEntry <- tkentry(inputFrame, text=tclVar(thin))
  thinLabel <- tklabel(inputFrame, text="thin")

  #-----------------------------------------------------------------------------
  # define buttons
  #-----------------------------------------------------------------------------
  runButton <- ttkbutton(lButtonFrame, width=12, text="Run", command=onRun)
  resetButton <- ttkbutton(lButtonFrame, width=12, text="Reset", command=onReset)
  cancelButton <- ttkbutton(lButtonFrame, width=12, text="Cancel", command=onCancel)

  nextplotButton <- ttkbutton(rButtonFrame, width=12, text="Next Plot", command=onNextplot)
  convButton <- ttkbutton(rButtonFrame, width=12, text="Converge", command=onConv)
  notconvButton <- ttkbutton(rButtonFrame, width=12, text="Not Converge", command=onNotconv)

  #-----------------------------------------------------------------------------
  ## tkgrid() and tkpack() the inputs and frames together
  #-----------------------------------------------------------------------------
  tkgrid(prior.lambda.lLabel, prior.lambda.l, prior.lambda.r, sticky="nw", padx=c(10,0), pady=c(0,15), columnspan=2)
  tkgrid(prior.pi.lLabel, prior.pi.l, prior.pi.r, sticky="nw", padx=c(10,0), pady=c(0,15), columnspan=2)
  tkgrid(chainsLabel, chainsEntry, sticky="nw", padx=c(10,0), pady=c(0,15), columnspan=3)
  tkgrid(burnLabel, burnEntry, sticky="nw", padx=c(10,0), pady=c(0,15), columnspan=3)
  tkgrid(updateLabel, updateEntry, sticky="nw", padx=c(10,0), pady=c(0,15), columnspan=3)
  tkgrid(thinLabel, thinEntry, sticky="nw", padx=c(10,0), pady=c(0,15), columnspan=3)
  tkgrid(runButton, resetButton, cancelButton, sticky="we", padx=c(5,5))
  tkgrid(nextplotButton, convButton, notconvButton, sticky="swe", padx=c(5,5))

  tkpack(inputFrame, side="top")
  tkpack(rButtonFrame, pady=c(5,25), side="bottom")
  tkpack(lButtonFrame, side="bottom", pady=c(0,25))
  tkpack(leftFrame, side="left")

  tkpack(imgFrame, side="top", padx=c(15,0), pady=c(0,10))
  tkpack(rightFrame, side="right")

  tkraise(bayesZIPWindow)

  tkwait.window(bayesZIPWindow) # otherwise, mod@convergence won't be saved

  #-----------------------------------------------------------------------------
  # prepare output
  #-----------------------------------------------------------------------------
  if(exists("mod", where=envirZIP)){
    mod <- get("mod", envir=envirZIP)
  } # end if
  return(mod)
} # end of function ZIPGUI



################################################################################
################################################################################
#' @details The diagnostic parameters \code{se} and \code{sp} are defined at
#' the pool level, equivalent to \code{missclass='pool'} in \code{rrisk.BayesPEM}
#' function. See \code{\link{rrisk.BayesPEM}} for more details.
#'
#' @description This function provides a GUI for the function \link{rrisk.BayesPEM}.
#'
#' @name PEMGUI
#' @aliases PEMGUI
#' @title GUI for Bayesian Prevalence estimation under misclassification (PEM)
#' @usage PEMGUI(x=20, n=20, k=10, prior.pi=c(1,19), prior.se=c(1,1),
#'  prior.sp=c(1,1), chains=3, burn=1000, update=10000, thin=1)
#' @param x scalar value for number of pools (\code{k>1}) or single outcomes (\code{k=1}) with positive test result
#' @param n scalar value for number of pools tested (\code{k>1}) or the sample size in application study (\code{k=1})
#' @param k scalar value for number of individual samples physically combined into one pool
#' @param prior.pi numeric vector containing parameters of a beta distribution as prior for prevalence \code{pi}, e.g. \code{pi} \eqn{\sim} \code{prior.pi(*,*)=beta(*,*)}
#' @param prior.se numeric vector containing parameters of a beta distribution as prior for sensitivity \code{se}, e.g. \code{se} \eqn{\sim} \code{prior.se(*,*)=beta(*,*)}
#' @param prior.sp numeric vector containing parameters of a beta distribution as prior for specificity \code{sp}, e.g. \code{sp} \eqn{\sim} \code{prior.sp(*,*)=beta(*,*)}
#' @param chains positive single numeric value, number of independent MCMC chains (default 3)
#' @param burn positive single numeric value, length of the burn-in period (default 1000)
#' @param update positive single numeric value, length of update iterations for estimation (default 10000)
#' @param thin positive single numeric value (default 1). The samples from every kth iteration will be used for
#'        inference, where k is the value of thin. Setting \code{thin > 1} can help to reduce the autocorrelation
#'        in the sample.
#' @return The function \code{PEMGUI} returns an instance of the \code{\linkS4class{bayesmodelClass}}
#' class containing following information
#' \item{\code{convergence}}{logical, whether the model has converged (assessed by the user)}
#' \item{\code{results}}{data frame containing statitsics of the posterior distribution}
#' \item{\code{jointpost}}{data frame giving the joint posterior probability distribution}
#' \item{\code{nodes}}{names of the parameters jointly estimated by the Bayes model}
#' \item{\code{model}}{model in BRugs/Winbugs syntax as a character string}
#' \item{\code{chains}}{number of independent MCMC chains}
#' \item{\code{burn}}{length of burn-in period}
#' \item{\code{update}}{length of update iterations for estimation}
#' @note The convergence of the model is assessed by the user using diagnostic plots
#' provided by the \pkg{BRugs} package.
#' @seealso \code{\link{rrisk.BayesPEM}}
#' @keywords manip
#' @export
# @references Greiner, M., Belgoroski, N., NN (2011). Estimating prevalence using diagnostic data
# with pooled samples: R function \code{rrisk.BayesPEM}. J.Stat.Software (in preparation).
#' @examples
#' #------------------------------------------
#' # Example of PEM model. Without parameters,
#' # the input fields will show default values
#' #------------------------------------------
#' \donttest{mod <- PEMGUI()}

PEMGUI <- function(x=20, n=20, k=10, prior.pi=c(1,19), prior.se=c(1,1),
                   prior.sp=c(1,1), chains=3,burn=1000, update=10000, thin=1){
  #-----------------------------------------------------------------------------
  # GUIDiag: function for plotting the two diagnosis plots
  #-----------------------------------------------------------------------------
  GUIDiag<-function(nodes,plotnumber=1)
  { len<-length(nodes)
  if(plotnumber==1){ par(mfrow=c(2,len))
    for(i in 1:len){
      if(len<=3){
        maintext<-paste("Gelman-Rubin conv. \nstatistic for",nodes[i])
      } else maintext<-paste("G.-R. statistic for",nodes[i])
      rrisk.samplesBgr(nodes[i],main=maintext,cex.main=1)
      abline(h=1,lty="dashed")
    }
    for(i in 1:len){
      plotDensity(nodes[i],main=paste("Density plot of ",nodes[i]),cex.main=1)
    }
  } else if (plotnumber==2){
    if(len<=3){
      par(mfrow=c(len,1))
    } else  par(mfrow=c(2,3))
    for(i in 1:len){
      plotHistory(nodes[i],main=paste("Convergence monitored for node",nodes[i]),cex.main=1)
    }
  } else if(plotnumber==3){
    if(len<=3){
      par(mfrow=c(len,1))
    } else  par(mfrow=c(2,3))
    for(i in 1:len){
      plotAutoC(nodes[i],main=paste("Autocorrelation plot for node",nodes[i]),cex.main=1)
    }
  }
  }
  #-----------------------------------------------------------------------------
  # what happens by pressing RUN button
  #-----------------------------------------------------------------------------
  onRun <- function(...){
    x. <- as.numeric(tclvalue(tkget(xEntry)))
    n. <- as.numeric(tclvalue(tkget(nEntry)))
    k. <- as.numeric(tclvalue(tkget(kEntry)))
    prior.pi.l. <- as.numeric(tclvalue(tkget(prior.pi.l)))
    prior.pi.r. <- as.numeric(tclvalue(tkget(prior.pi.r)))
    prior.se.l. <- as.numeric(tclvalue(tkget(prior.se.l)))
    prior.se.r. <- as.numeric(tclvalue(tkget(prior.se.r)))
    prior.sp.l. <- as.numeric(tclvalue(tkget(prior.sp.l)))
    prior.sp.r. <- as.numeric(tclvalue(tkget(prior.sp.r)))
    chains. <- as.numeric(tclvalue(tkget(chainsEntry)))
    burn. <- as.numeric(tclvalue(tkget(burnEntry)))
    update. <-  as.numeric(tclvalue(tkget(updateEntry)))
    thin. <-  as.numeric(tclvalue(tkget(thinEntry)))
    #misclass. <- misclasses[as.numeric(tclvalue(tcl(misclassEntry, "getvalue")))+1]

    # testmodel:
    # rrisk.BayesPEM(20, 20, 10, simulation=FALSE, prior.pi=c(0.1,0.2), prior.se=c(0.1,0.2), prior.sp=c(0.1,0.2), chains=10, misclass="pool")

    mod <- rrisk.BayesPEM(x=x., n=n., k=k., simulation=TRUE, prior.pi=c(prior.pi.l., prior.pi.r.),
                          prior.se=c(prior.se.l., prior.se.r.), prior.sp=c(prior.sp.l., prior.sp.r.),
                          chains=chains., misclass="pool",update=update., burn=burn.,thin=thin.)

    assign("mod", value=mod, envir=envirPEM)

    if(exists("imgPlot1", where=envirPEM)){
      imgPlot1 <- get("imgPlot1", envir=envirPEM)
      tkpack.forget(imgPlot1)
    }
    if(exists("imgPlot2", where=envirPEM)){
      imgPlot2 <- get("imgPlot2", envir=envirPEM)
      tkpack.forget(imgPlot2)
    }
    if(exists("imgPlot3", where=envirPEM)){
      imgPlot3 <- get("imgPlot3", envir=envirPEM)
      tkpack.forget(imgPlot3)
    }

    nodes <- c("pi", "se", "sp")
    imgPlot1 <- tkrplot(imgFrame,fun=function()GUIDiag(nodes,plotnumber=1),hscale=1.65,vscale=1.5)
    imgPlot2 <- tkrplot(imgFrame,fun=function()GUIDiag(nodes,plotnumber=2),hscale=1.55,vscale=1.7)
    imgPlot3 <- tkrplot(imgFrame,fun=function()GUIDiag(nodes,plotnumber=3),hscale=1.55,vscale=1.7)

    assign("imgPlot1", value=imgPlot1, envir=envirPEM)
    assign("imgPlot2", value=imgPlot2, envir=envirPEM)
    assign("imgPlot3", value=imgPlot3, envir=envirPEM)

    tkpack(imgPlot1, side="top")
    assign("plotNumber",value=1,envir=envirPEM)

    tkraise(bayesPEMWindow)
  }

  #-----------------------------------------------------------------------------
  # what happens by pressing RESET button
  #-----------------------------------------------------------------------------
  onReset <- function(...){
    tkconfigure(xEntry, text=tclVar(x))
    tkconfigure(nEntry, text=tclVar(n))
    tkconfigure(kEntry, text=tclVar(k))
    tkconfigure(prior.pi.l, text=tclVar(prior.pi[1]))
    tkconfigure(prior.pi.r, text=tclVar(prior.pi[2]))
    tkconfigure(prior.se.l, text=tclVar(prior.se[1]))
    tkconfigure(prior.se.r, text=tclVar(prior.se[2]))
    tkconfigure(prior.sp.l, text=tclVar(prior.sp[1]))
    tkconfigure(prior.sp.r, text=tclVar(prior.sp[2]))
    #tkconfigure(misclassEntry, text=misclass)
    tkconfigure(chainsEntry, text=tclVar(chains))
    tkconfigure(burnEntry, text=tclVar(burn))
    tkconfigure(updateEntry, text=tclVar(update))
    tkconfigure(thinEntry, text=tclVar(thin))
  }

  #-----------------------------------------------------------------------------
  # what happens by pressing CANCEL button
  #-----------------------------------------------------------------------------
  onCancel <- function(...){
    tkdestroy(bayesPEMWindow)
  }

  #-----------------------------------------------------------------------------
  # what happens by pressing NEXT button
  #-----------------------------------------------------------------------------
  onNextplot <- function(...){
    # get the current plot number (in 1,2) and switch to the other plot
    if(!exists("plotNumber", envirPEM)) return(FALSE) # if there was no model run yet

    plotNumber <- get("plotNumber", envir=envirPEM)
    if(plotNumber==1){
      imgPlot1 <- get("imgPlot1", envir=envirPEM)
      tkpack.forget(imgPlot1)
      assign("plotNumber",value=2,envir=envirPEM)
      imgPlot2 <- get("imgPlot2", envir=envirPEM)
      tkpack(imgPlot2, side="top")
    } else if(plotNumber==2){
      imgPlot2 <- get("imgPlot2", envir=envirPEM)
      tkpack.forget(imgPlot2)
      assign("plotNumber",value=3,envir=envirPEM)
      imgPlot3 <- get("imgPlot3", envir=envirPEM)
      tkpack(imgPlot3, side="top")
    } else if(plotNumber==3){
      imgPlot3 <- get("imgPlot3", envir=envirPEM)
      tkpack.forget(imgPlot3)
      assign("plotNumber",value=1,envir=envirPEM)
      imgPlot1 <- get("imgPlot1", envir=envirPEM)
      tkpack(imgPlot1, side="top")
    }
    tkraise(bayesPEMWindow)
  }

  #-----------------------------------------------------------------------------
  # what happens by pressing CONVERGE button
  #-----------------------------------------------------------------------------
  onConv <- function(...){
    mod <- get("mod", envir=envirPEM)
    mod@convergence <- TRUE
    assign("mod", value=mod, envir=envirPEM)
    tkdestroy(bayesPEMWindow)
  }

  #-----------------------------------------------------------------------------
  # what happens by pressing NOT CONVERGE button
  #-----------------------------------------------------------------------------
  onNotconv <- function(...){
    mod <- get("mod", envir=envirPEM)
    mod@convergence <- FALSE
    assign("mod", value=mod, envir=envirPEM)
    tkdestroy(bayesPEMWindow)
  }

  #-----------------------------------------------------------------------------
  # define Dialof window
  #-----------------------------------------------------------------------------
  tclRequire("BWidget")

  assign("envirPEM",value=new.env(),envir=.GlobalEnv)

  mod <- new("bayesmodelClass")

  #-----------------------------------------------------------------------------
  ## define GUI window and frames
  #-----------------------------------------------------------------------------
  bayesPEMWindow <- tktoplevel()
  tkwm.title(bayesPEMWindow, "GUI for the function rrisk.bayesPEM")
  tkwm.resizable(bayesPEMWindow, 0, 0) # fixed size
  tkwm.maxsize(bayesPEMWindow,900,600)
  tkwm.minsize(bayesPEMWindow,900,600)

  leftFrame <- tkframe(bayesPEMWindow)
  rightFrame <- tkframe(bayesPEMWindow)
  imgFrame <- tkframe(rightFrame, height=600,width=200)
  inputFrame <- tkframe(leftFrame)
  lButtonFrame <- tkframe(leftFrame)
  rButtonFrame <- tkframe(leftFrame)

  #-----------------------------------------------------------------------------
  ## define input fields
  #-----------------------------------------------------------------------------
  xEntry <- tkentry(inputFrame, text=tclVar(x))
  xLabel <- tklabel(inputFrame, text="x")

  nEntry <- tkentry(inputFrame, text=tclVar(n))
  nLabel <- tklabel(inputFrame, text="n")

  kEntry <- tkentry(inputFrame, text=tclVar(k))
  kLabel <- tklabel(inputFrame, text="k")

  prior.pi.l <- tkentry(inputFrame, text=tclVar(prior.pi[1]), width=6)
  prior.pi.lLabel <- tklabel(inputFrame, text="prior.pi (beta)")
  prior.pi.r <- tkentry(inputFrame, text=tclVar(prior.pi[2]), width=6)

  prior.se.l <- tkentry(inputFrame, text=tclVar(prior.se[1]), width=6)
  prior.se.lLabel <- tklabel(inputFrame, text="prior.se (beta)")
  prior.se.r <- tkentry(inputFrame, text=tclVar(prior.se[2]), width=6)

  prior.sp.l <- tkentry(inputFrame, text=tclVar(prior.sp[1]), width=6)
  prior.sp.lLabel <- tklabel(inputFrame, text="prior.sp (beta)")
  prior.sp.r <- tkentry(inputFrame, text=tclVar(prior.sp[2]), width=6)

  #misclassLabel <- tklabel(inputFrame, text="misclass")
  #misclasses <- c("pool", "individual", "compare")
  #misclassEntry <- tkwidget(inputFrame, "ComboBox", editable=FALSE, values=misclasses, text=misclass)

  chainsEntry <- tkentry(inputFrame, text=tclVar(chains))
  chainsLabel <- tklabel(inputFrame, text="chains")

  burnEntry <- tkentry(inputFrame, text=tclVar(burn))
  burnLabel <- tklabel(inputFrame, text="burn")

  updateEntry <- tkentry(inputFrame, text=tclVar(update))
  updateLabel <- tklabel(inputFrame, text="update")

  thinEntry <- tkentry(inputFrame, text=tclVar(thin))
  thinLabel <- tklabel(inputFrame, text="thin")

  runButton <- ttkbutton(lButtonFrame, width=12, text="Run", command=onRun)
  resetButton <- ttkbutton(lButtonFrame, width=12, text="Reset", command=onReset)
  cancelButton <- ttkbutton(lButtonFrame, width=12, text="Cancel", command=onCancel)

  nextplotButton <- ttkbutton(rButtonFrame, width=12, text="Next Plot", command=onNextplot)
  convButton <- ttkbutton(rButtonFrame, width=12, text="Converge", command=onConv)
  notconvButton <- ttkbutton(rButtonFrame, width=12, text="Not Converge", command=onNotconv)


  ## tkgrid() and tkpack() the inputs and frames together
  ## must be run "inside-out"
  ################################################################

  tkgrid(xLabel, xEntry, sticky="nw", padx=c(10,10), pady=c(10,15), columnspan=3)
  tkgrid(nLabel, nEntry, sticky="nw", padx=c(10,10), pady=c(0,15), columnspan=3)
  tkgrid(kLabel, kEntry, sticky="nw", padx=c(10,10), pady=c(0,15), columnspan=3)
  tkgrid(prior.pi.lLabel, prior.pi.l, prior.pi.r, sticky="nw", padx=c(10,0), pady=c(0,15), columnspan=2)
  tkgrid(prior.se.lLabel, prior.se.l, prior.se.r, sticky="nw", padx=c(10,0), pady=c(0,15), columnspan=2)
  tkgrid(prior.sp.lLabel, prior.sp.l, prior.sp.r, sticky="nw", padx=c(10,0), pady=c(0,15), columnspan=2)
  #tkgrid(misclassLabel, misclassEntry, columnspan=3, padx=c(10,0), pady=c(0,15))
  tkgrid(chainsLabel, chainsEntry, sticky="nw", padx=c(10,0), pady=c(0,15), columnspan=3)
  tkgrid(burnLabel, burnEntry, sticky="nw", padx=c(10,0), pady=c(0,15), columnspan=3)
  tkgrid(updateLabel, updateEntry, sticky="nw", padx=c(10,0), pady=c(0,15), columnspan=3)
  tkgrid(thinLabel, thinEntry, sticky="nw", padx=c(10,0), pady=c(0,15), columnspan=3)

  tkgrid(runButton, resetButton, cancelButton, sticky="we", padx=c(5,5))
  tkgrid(nextplotButton, convButton, notconvButton, sticky="swe", padx=c(5,5))

  tkpack(inputFrame, side="top")
  tkpack(rButtonFrame, pady=c(5,25), side="bottom")
  tkpack(lButtonFrame, side="bottom", pady=c(0,25))
  tkpack(leftFrame, side="left")

  tkpack(imgFrame, side="top", padx=c(15,0), pady=c(0,10))
  tkpack(rightFrame, side="right") # zum zweiten mal

  tkraise(bayesPEMWindow)

  tkwait.window(bayesPEMWindow) # otherwise, mod@convergence won't be saved

  if(exists("mod", where=envirPEM)){
    mod <- get("mod", envir=envirPEM)
  } # end if statement

  return(mod)
} # end of function PEMGUI

