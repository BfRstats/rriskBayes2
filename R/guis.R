
# pressing the next button ------------------------------------------------

onNextplot <- function(envir, nvars, ...){
  # get the current plot number (in 1,2) and switch to the other plot
  if(!exists("plotNumber", envir)) return(FALSE) # if there was no model run yet
  
  plotNumber <- get("plotNumber", envir=envir)
  
  if(plotNumber <= nvars){
    imgPlot1 <- get(paste0("imgPlot", plotNumber), envir=envir)
    tkpack.forget(imgPlot1)
    assign("plotNumber",value=plotNumber+1,envir=envir)
    imgPlot2 <- get(paste0("imgPlot",plotNumber+1), envir=envir)
    tkpack(imgPlot2, side="top")
  }else if(plotNumber == nvars + 1){       
    imgPlot3 <- get(paste0("imgPlot", plotNumber), envir=envir)
    tkpack.forget(imgPlot3)
    assign("plotNumber",value=1,envir=envir)
    imgPlot1 <- get("imgPlot1", envir=envir)
    tkpack(imgPlot1, side="top")
  }
} 


# ZIPGUI ------------------------------------------------------------------

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


ZIPGUI <- function(data, prior.lambda=c(0, 100), prior.pi=c(1, 1),
                   chains=3, burn=1000, update=10000,thin=1){
  if(missing(data))
    stop("please provide data of minimal length 10!")
  
  # what happens by pressing RUN button
  
  onRun <- function(...){
    prior.pi.l. <- as.numeric(tclvalue(tkget(prior.pi.l)))
    prior.pi.r. <- as.numeric(tclvalue(tkget(prior.pi.r)))
    prior.lambda.l. <- as.numeric(tclvalue(tkget(prior.lambda.l)))
    prior.lambda.r. <- as.numeric(tclvalue(tkget(prior.lambda.r)))
    chains. <- as.numeric(tclvalue(tkget(chainsEntry)))
    burn. <- as.numeric(tclvalue(tkget(burnEntry)))
    update. <- as.numeric(tclvalue(tkget(updateEntry)))
    thin. <- as.numeric(tclvalue(tkget(thinEntry)))
    
    mod <- try(rrisk.BayesZIP(data=data, prior.pi=c(prior.pi.l., prior.pi.r.), 
                              prior.lambda=c(prior.lambda.l., prior.lambda.r.),
                              chains=chains., burn=burn., update=update.,thin=thin.))
    
    if(inherits(mod, "try-error")){
      tkdestroy(bayesZIPWindow) 
      stop("Any error occured during fitting bayesian ZIP model",call.=FALSE)
    }
    assign("mod", value=mod, envir=envirZIP)
    vars <- mod@results$monitor
    nvars <- length(vars)
    assign("nvars", nvars, envir=envirZIP)
 
    #"trace", "ecdf", "histogram", "autocorr" plots
    k <- 1
    for(i in 1:nvars) {
      assign(paste0("imgPlot", k), value = tkrplot(imgFrame,fun=function() plot(mod@results, vars=vars[k]), hscale=1.6, vscale=1.5), envir=envirZIP)
      k <- k+1
    }
    #"crosscorrelation" plot
    assign(paste0("imgPlot", nvars+1), value = tkrplot(imgFrame,fun=function() plot(mod@results, vars=vars, plot.type="crosscorr"), hscale=1.6, vscale=1.5), envir=envirZIP)
 
    imgPlot1 <- get("imgPlot1", envir=envirZIP)
    tkpack(imgPlot1, side="top")
    
    assign("plotNumber", value=1, envir = envirZIP)
    
    tkraise(bayesZIPWindow)
    
  } # end onRun() fucntion
  
  
  # what happens by pressing RESET button
  
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
  
  
  # what happens by pressing CANCEL button
  
  onCancel <- function(...){
    tkdestroy(bayesZIPWindow) 
  } # end onCancel() function
  
  
  # what happens by pressing Nexplot button
  
  
  # what happens by pressing CONVERGENCE button
  
  onConv <- function(...){
    mod <- get("mod", envir=envirZIP)
    mod@convergence <- TRUE
    assign("mod", value=mod, envir=envirZIP)
    tkdestroy(bayesZIPWindow)
  } # end fucntion onConv()
  
  
  # what happens by pressing CONVERGENCE button
  
  onNotconv <- function(...){
    mod <- get("mod", envir=envirZIP)
    mod@convergence <- FALSE
    assign("mod", value=mod, envir=envirZIP)
    tkdestroy(bayesZIPWindow)
  } # end fucntion onNotconv()
  
  
  # define help varriable(s)
  
  assign("envirZIP", value=new.env(), envir=.GlobalEnv)
  
  
  # define GUI window and frames
  
  bayesZIPWindow <- tktoplevel()
  tkwm.title(bayesZIPWindow, "GUI for the function rrisk.bayesZIP")
  tkwm.resizable(bayesZIPWindow, 0, 0) # fixed size
  tkwm.maxsize(bayesZIPWindow,1000,600)
  tkwm.minsize(bayesZIPWindow,1000,600 )
  
  mod <- new("bayesmodelClass")
  
  leftFrame <- tkframe(bayesZIPWindow)
  rightFrame <- tkframe(bayesZIPWindow)
  imgFrame <- tkframe(rightFrame,height=600,width=200)
  inputFrame <- tkframe(leftFrame)
  lButtonFrame <- tkframe(leftFrame)
  rButtonFrame <- tkframe(leftFrame)
  #rButtonFrame <- tkframe(rightFrame)
  
  
  # define input fields
  
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
  
  
  # define buttons
  
  runButton <- ttkbutton(lButtonFrame, width=12, text="Run", command=onRun)
  resetButton <- ttkbutton(lButtonFrame, width=12, text="Reset", command=onReset)
  cancelButton <- ttkbutton(lButtonFrame, width=12, text="Cancel", command=onCancel)

  nextplotButton <- ttkbutton(rButtonFrame, width=12, text="Next Plot", command=function(...) onNextplot(envir=envirZIP, nvars=2))
  convButton <- ttkbutton(rButtonFrame, width=12, text="Converge", command=onConv)
  notconvButton <- ttkbutton(rButtonFrame, width=12, text="Not Converge", command=onNotconv)
  
  
  ## tkgrid() and tkpack() the inputs and frames together
  
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
  
  
  # prepare output
  
  if(exists("mod", where=envirZIP)){
    mod <- get("mod", envir=envirZIP) 
  } # end if
  return(mod)
} # end of function ZIPGUI



# PEMGUI ------------------------------------------------------------------

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

PEMGUI <- function(x=20, n=20, k=1, prior.pi=c(1,1), prior.se=c(1,1),
                   prior.sp=c(1,1), chains=3, burn=1000, update=10000, thin=1){
  
  # GUIDiag: function for plotting the two diagnosis plots
  
  
  
  # what happens by pressing RUN button
  
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
    
    mod <- try(rrisk.BayesPEM(x=x., n=n., k=k., prior.pi=c(prior.pi.l., prior.pi.r.), 
                              prior.se=c(prior.se.l., prior.se.r.), prior.sp=c(prior.sp.l., prior.sp.r.),
                              chains=chains., misclass="pool",update=update., burn=burn.,thin=thin.))
    
    if(inherits(mod, "try-error")){
      tkdestroy(bayesPEMWindow)
      stop("Any error occured during fitting bayesian PEM model",call.=FALSE)
    } 
    
    assign("mod", value=mod, envir=envirPEM)
    
    vars <- mod@results$monitor
    nvars <- length(vars)
    assign("nvars", nvars, envir=envirPEM)
   
    #"trace", "ecdf", "histogram", "autocorr" plots
    k <- 1
    for(i in 1:nvars) {
      assign(paste0("imgPlot", k), value = tkrplot(imgFrame,fun=function() plot(mod@results, vars=vars[k]), hscale=1.6, vscale=1.5), envir=envirPEM)
      k <- k+1
    }
    #crosscorrelation plot
    assign(paste0("imgPlot", nvars+1), value = tkrplot(imgFrame,fun=function() plot(mod@results, vars=vars, plot.type="crosscorr"), hscale=1.6, vscale=1.5), envir=envirPEM)
    
    imgPlot1 <- get("imgPlot1", envir=envirPEM)
    tkpack(imgPlot1, side="top")
    assign("plotNumber",value=1,envir=envirPEM)
    tkraise(bayesPEMWindow)
  }
  
  
  # what happens by pressing RESET button
  
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
  
  
  # what happens by pressing CANCEL button
  
  onCancel <- function(...){
    tkdestroy(bayesPEMWindow)
  }
  
  
  # what happens by pressing Nexplot button - see common function 
  

  
  # what happens by pressing CONVERGE button
  
  onConv <- function(...){
    mod <- get("mod", envir=envirPEM)
    mod@convergence <- TRUE
    assign("mod", value=mod, envir=envirPEM)
    tkdestroy(bayesPEMWindow)
  }
  
  
  # what happens by pressing NOT CONVERGE button
  
  onNotconv <- function(...){
    mod <- get("mod", envir=envirPEM)
    mod@convergence <- FALSE
    assign("mod", value=mod, envir=envirPEM)
    tkdestroy(bayesPEMWindow)
  }
  
  
  # define Dialof window
   
  tclRequire("BWidget")
  
  assign("envirPEM",value=new.env(),envir=.GlobalEnv)
  
  mod <- new("bayesmodelClass")
  
  
  ## define GUI window and frames
  
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
  
  
  ## define input fields
  
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

  nextplotButton <- ttkbutton(rButtonFrame, width=12, text="Next Plot", command=function(...) onNextplot(envir=envirPEM, nvars=3))
  convButton <- ttkbutton(rButtonFrame, width=12, text="Converge", command=onConv)
  notconvButton <- ttkbutton(rButtonFrame, width=12, text="Not Converge", command=onNotconv)
  
  
  ## tkgrid() and tkpack() the inputs and frames together
  ## must be run "inside-out"
  
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
  } 
  return(mod)
} # end of function PEMGUI


