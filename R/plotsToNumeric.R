##from runjags:::runjagsplots in runjags:::plot.runjags 

#res <- test_rrisk.BayesPEM("individual")

plots2num <- function(){
  
  separate.chains <- TRUE
  mcmclist <- res@results$mcmc  
  thinned.mcmc <- combine.mcmc(list(mcmclist), collapse.chains = FALSE, 
                               return.samples = niter(mcmclist))
  nchains <- nchain(thinned.mcmc)
  chainnames <- paste("chain_", gsub(" ", "0", 
                                     format(1:nchains, scientific = FALSE)), sep = "")
  niters <- niter(thinned.mcmc)
  plot1 = plot2 = plot3 = plot4 = plot5 <- lapply(1:length(varnames(thinned.mcmc)), 
                                                  function(x) {
                                                    newl <- vector("list", length = nchains)
                                                    names(newl) <- chainnames
                                                    return(newl)
                                                  })
  names(plot1) = names(plot2) = names(plot3) = names(plot4) = names(plot5) <- varnames(thinned.mcmc)
  dataAutocorr <-  dataECDF <- dataDens <-  dataTrace <- list()
  
  
  for (i in 1:length(varnames(thinned.mcmc))) {
    all.thinned <- as.mcmc.list(lapply(thinned.mcmc, 
                                       function(x) return(x[, i, drop = FALSE])))
    all.full <- as.mcmc.list(lapply(mcmclist, function(x) return(x[,i, drop = FALSE])))
    varname <- dimnames(all.thinned[[1]])[[2]]
    
    dataAutocorr[[i]] <-  dataECDF[[i]] <- dataDens[[i]] <-  dataTrace[[i]] <- list()
    
    
    for (p in 1:nchains) {
      thinplotdata <- all.thinned[[p]]
      allplotdata <- all.full[[p]]
      chain <- p
     # nchains <- 1
      
      ##trace###################################################################
      plotopts <- list()
      dataTrace[[i]][[p]] <- thinplotdata
      plotopts$x <- thinplotdata
      plot1[[i]][[p]] <- do.call("xyplot", args = plotopts)
      
      ##density################################################################
      plotopts <- list()
      dataDens[[i]][[p]] <- allplotdata
      plotopts$x <- allplotdata
      plot2[[i]][[p]] <- do.call("densityplot", 
                                 args = plotopts)
      
      ##ecdf####################################################################
      plotopts <- list()
      tdat <- data.frame(cd = 1:niters/niters, 
                         vx = numeric(niters * max(1, nchains * 
                                                     (!separate.chains))))
      tdat$vx <- quantile(as.numeric(unlist(allplotdata)), 
                          probs = tdat$cd)
      plotopts$data <- tdat
      plotopts$x <- cd ~ vx
      
      dataECDF[[i]][[p]] <- tdat
      plot3[[i]][[p]] <- do.call("xyplot", args = plotopts)
      
      plotdata <- allplotdata
      
      ##autocorr################################################################
      plotopts <- list()
      acfs <- acf(plotdata, plot = FALSE)
      plotopts$data <- data.frame(acfs = as.numeric(acfs$acf), 
                                  lags = as.numeric(acfs$lag))
      dataAutocorr[[i]][[p]] <- plotopts$data
      plotopts$x <- acfs ~ lags
      plot5[[i]][[p]] <- do.call("barchart", args = plotopts)
      
    }
  }
  ##crosscorr
  ccplot <- list(crosscorr = crosscorr.plot(mcmclist))
  class(ccplot) <- "runjagsplots"
  
  return(list(dataAutocorr=dataAutocorr, dataECDF=dataECDF, dataDens=dataDens, dataTrace=dataTrace))
}

#kruskal.test(resnum$dataTrace[[1]])
#a <- acf(resnum$dataTrace[[1]][[1]], lag=4)


##pairwise kolmogorov smirnov
# data <- do.call(cbind, resnum$dataTrace[[1]])
# pwKS <- function(x, y, ..., 
#               alternative = c("two.sided"), exact = NULL){ 
#   #w <- getOption("warn") 
#   #options(warn = -1)  # ignore warnings 
#   p <- ks.test(x, y, ..., alternative = alternative, exact = 
#                  exact)$p.value 
#   #options(warn = w) 
#   p 
# } 
# pvals <- apply(data, 2, function(x) apply(data, 2, function(y) pwKS(x, y)))
# p <- as.numeric(pvals)
# p.adjust(p=p, method = p.adjust.methods, n = length(p))

