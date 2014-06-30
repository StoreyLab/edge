#' Fitted models for an \code{\linkS4class{edgeFit}} object.
#'
#' @param object \code{\linkS4class{edgeFit}} object
#' 
#' @usage modelFits(object)
#' 
#' @return \code{modelFits} returns the formula for the fitted models.
#'
#' @author John Storey, Andrew Bass
#' 
#' @seealso \code{\link{edgeFit}}, \code{\link{edgeSet-class}}
#' 
#' @keywords modelFits
#' 
#' @rdname modelFits-method
#' 
#' @exportMethod modelFits
setGeneric("modelFits", function(object) standardGeneric("modelFits"))

#' Regression coefficients from full model fit
#'
#' Access the full model fitted coefficients of an \code{\linkS4class{edgeFit}} object.
#'
#' @param object \code{\linkS4class{edgeFit}}
#' 
#' @usage betaCoef(object)
#' 
#' @return \code{betaCoef} returns the regression coefficients.
#'
#' @author John Storey, Andrew Bass 
#' 
#' @seealso \code{\link{edgeFit}}
#' 
#' @keywords betaCoef
#' 
#' @rdname betaCoef-method
#' 
#' @exportMethod betaCoef
setGeneric("betaCoef", function(object) standardGeneric("betaCoef"))

#' Statistical method used in analysis
#'
#' Access the statistic type in an \code{\linkS4class{edgeFit}} object. Can either be the Optimal Discovery Procedure (odp) or the likelihood ratio test (lrt).
#'
#' @param object \code{\linkS4class{edgeFit}}
#' 
#' @usage sType(object)
#' 
#' @return \code{sType} returns the statistic type- either "odp" or "lrt".
#'
#' @author John Storey, Andrew Bass 
#' 
#' @seealso \code{\link{edgeFit}}, \code{\linkS4class{edgeFit}} and 
#' \code{\linkS4class{edgeSet}}
#' 
#' @keywords sType
#' 
#' @rdname sType-method
#' 
#' @exportMethod sType
setGeneric("sType", function(object) standardGeneric("sType"))

#' Fitted data from the full model 
#'
#' Access the fitted data from the full model in an \code{\linkS4class{edgeFit}} object.
#'
#' @param object \code{\linkS4class{edgeFit}}
#'  
#' @usage fitFull(object, ...)
#' 
#' @return \code{fitFull} returns a matrix of fitted values from full model.
#'
#' @author John Storey, Andrew Bass
#' 
#' @seealso \code{\link{edgeFit}}
#' 
#' @keywords fitFull
#' 
#' @rdname fitFull-method
#' 
#' @exportMethod fitFull
setGeneric("fitFull", function(object, ...) standardGeneric("fitFull"))

#' Fitted data from the null model
#'
#' Access the fitted data from the null model in an \code{\linkS4class{edgeFit}} object.
#'
#' @param object \code{\linkS4class{edgeFit}}
#' 
#' @usage fitNull(object, ...)
#' 
#' @return \code{fitFull} returns a matrix of fitted values from null model.
#'
#' @author  John Storey, Andrew Bass
#' 
#' @seealso \code{\link{edgeFit}}
#' 
#' @keywords fitNull
#' 
#' @rdname fitNull-method
#' 
#' @exportMethod fitNull  
setGeneric("fitNull", function(object, ...) standardGeneric("fitNull"))

#' Residuals of full model fit
#'
#' Access the fitted full model residuals in an \code{\linkS4class{edgeFit}} object.
#'
#' @param object \code{\linkS4class{edgeFit}}
#' 
#' @usage resFull(object)
#' 
#' @return \code{resFull} returns a matrix of residuals from full model.
#'
#' @author John Storey, Andrew Bass
#' 
#' @seealso \code{\link{edgeFit}}
#' 
#' @keywords resFull
#' 
#' @rdname resFull-method
#' 
#' @exportMethod resFull
setGeneric("resFull", function(object) standardGeneric("resFull"))

#' Residuals of null model fit
#'
#' Access the fitted null model residuals in an \code{\linkS4class{edgeFit}} object.
#'
#' @param object \code{\linkS4class{edgeFit}}
#' 
#' @usage resNull(object)
#' 
#' @return \code{resNull} returns a matrix of residuals from null model.
#'
#' @author John Storey, Andrew Bass
#' 
#' @seealso  \code{\link{edgeFit}}
#' 
#' @keywords  resNull
#' 
#' @rdname resNull-method
#' 
#' @exportMethod resNull 
setGeneric("resNull", function(object) standardGeneric("resNull"))
