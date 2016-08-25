# bayesmodelClass ---------------------------------------------------------


#' @description S4 class for displaying the output of the functions \code{rrisk.BayesPEM} and
#' \code{rrisk.BayesZIP}.
#'
#' @name bayesmodelClass-class
#' @aliases bayesmodelClass
#' @docType class
#' @title S4 class for displaying output of Bayesian models
#' @slot convergence {logical, whether the model has converged (assessed by the user only if GUI was used)}
#' @slot results summary of runjags object
#' @slot jointposterior {data frame giving the joint posterior probability distribution}
#' @slot nodes {names of the parameters jointly estimated by the Bayes model}
#' @slot model {model in rjags/JAGS syntax}
#' @slot chains {number of independent MCMC chains}
#' @slot burn {length of burn-in period}
#' @slot update {length of update iterations for estimation}
#' @rdname bayesmodelClass-class
#' @exportClass bayesmodelClass
#' @examples
#' new("bayesmodelClass")
#' new("bayesmodelClass",model="Some model...", nodes="Nodes info...")

setClass(Class="bayesmodelClass",
         representation=representation(
           convergence="logical",
           results="ANY",
           jointpost="ANY",
           nodes="ANY",
           model="ANY",
           chains="ANY",
           burn="ANY",
           update="ANY"),
         prototype=prototype(
           convergence=FALSE,
           results=NULL,
           jointpost=NULL,
           nodes=NULL,
           model=NULL,
           chains=NULL,
           burn=NULL,
           update=NULL))



# methods show  -----------------------------------------------------------

#' @description Show method for \code{\linkS4class{bayesmodelClass}}
#'
#' @name show-methods
#' @aliases show,bayesmodelClass-method
#' @docType methods
#' @title Show method for bayesmodelClass
#' @param object a \code{bayesmodelClass} object
#' @exportMethod show
#' @importFrom methods show
#' @examples
#' show(new("bayesmodelClass"))

setMethod(f="show",
          signature=signature(object="bayesmodelClass"),
          definition=function(object)
          { 
            cat("\n")
            if(is.null(object@convergence))
            { cat("convergence: \n NULL\n\n")
            } else cat("convergence: \n",object@convergence, "\n\n")
            if(is.null(object@results))
            { cat("results: \n NULL\n\n")
            } else {cat("results: \n");print(object@results);cat("\n")}
            if(is.null(object@jointpost))
            { cat("jointpost: \n NULL\n\n")
            } else {cat("jointpost: \n");print(head(object@jointpost));cat( "........\n\n")}
            if(is.null(object@nodes))
            { cat("nodes: \n NULL\n\n")
            } else cat("nodes: \n",object@nodes, "\n\n")
            if(is.null(object@model))
            { cat("model: \n NULL\n\n")
            } else cat("model: \n",object@model, "\n\n")
            if(is.null(object@chains))
            { cat("chains: \n NULL\n\n")
            } else cat("chains: \n",object@chains, "\n\n")
            if(is.null(object@burn))
            { cat("burn: \n NULL\n\n")
            } else cat("burn: \n",object@burn, "\n\n")
            if(is.null(object@update))
            { cat("update: \n NULL\n\n")
            } else cat("update: \n",object@update, "\n\n")})



