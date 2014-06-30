#' Full model equation from edgeSet object
#' 
#' These generic functions access and set the full model for 
#' \code{\linkS4class{edgeSet}} object.
#'
#' @param object \code{\linkS4class{edgeSet}}
#' @param value \code{formula}: The experiment design for the full model.
#' 
#' @usage fullModel(object)
#' 
#' @return \code{fullModel} returns the formula for the full model.
#'    
#' @author John Storey, Andrew Bass
#' 
#' @seealso \code{\linkS4class{edgeSet}}
#' 
#' @keywords fullModel, fullModel<-
#' 
#' @rdname fullModel-method
#' 
#' @exportMethod fullModel
setGeneric("fullModel", function(object) standardGeneric("fullModel"))

#' @rdname fullModel-method
#' 
#' @exportMethod fullModel<-
setGeneric("fullModel<-", function(object, value) {
  standardGeneric("fullModel<-") 
})

#' Null model equation from edgeSet object
#'
#' These generic functions access and set the null model for 
#' \code{\linkS4class{edgeSet}} object.
#'
#' @param object \code{\linkS4class{edgeSet}}
#' @param value \code{formula}: The experiment design for the null model.
#' 
#' @usage nullModel(object)
#' 
#' @return \code{nullModel} returns the formula for the null model.
#'
#' @author John Storey, Andrew Bass
#' 
#' @seealso \code{\linkS4class{edgeSet}}
#' 
#' @keywords nullModel, nullModel<-
#' 
#' @rdname nullModel-method
#' 
#' @exportMethod nullModel
setGeneric("nullModel", function(object) standardGeneric("nullModel"))

#' @rdname nullModel-method
#' 
#' @exportMethod nullModel<-
setGeneric("nullModel<-", function(object, value) {
  standardGeneric("nullModel<-") 
})

#' Matrix representation of null model
#'
#' These generic functions access and set the null matrix for 
#' \code{\linkS4class{edgeSet}} object.
#'
#' @param object \code{\linkS4class{edgeSet}}
#' @param value \code{matrix}: null model matrix where columns are covariates and rows are observations
#' 
#' @usage nullMatrix(object)
#' 
#' @return \code{nullMatrix} returns the value of the null model matrix.
#'
#' @author John Storey, Andrew Bass
#' 
#' @seealso \code{\linkS4class{edgeSet}}, \code{\link{fullModel}} and 
#' \code{\link{fullModel}}
#' 
#' @keywords nullMatrix, nullMatrix<-
#' 
#' @rdname nullMatrix-method   
#' 
#' @exportMethod nullMatrix
setGeneric("nullMatrix", function(object) standardGeneric("nullMatrix"))

#' @rdname nullMatrix-method
#'  
#' @exportMethod nullMatrix<-
setGeneric("nullMatrix<-", function(object, value) {
  standardGeneric("nullMatrix<-") 
})

#' Matrix representation of full model
#'
#' These generic functions access and set the full matrix for 
#' \code{\linkS4class{edgeSet}} object.
#'
#' @param object \code{\linkS4class{edgeSet}}
#' @param value \code{matrix}: full model matrix where the columns are the covariates amd rows are observations
#' 
#' @usage fullMatrix(object)
#' 
#' @return \code{fullMatrix} returns the value of the full model matrix.
#'
#' @author Andrew Bass, John Storey
#' 
#' @seealso \code{\linkS4class{edgeSet}}, \code{\link{fullModel}}
#' 
#' @keywords fullMatrix, fullMatrix<-
#' 
#' @rdname fullMatrix-method
#' 
#' @exportMethod fullMatrix
setGeneric("fullMatrix", function(object) standardGeneric("fullMatrix"))

#' @rdname fullMatrix-method
#' 
#' @exportMethod fullMatrix<-
setGeneric("fullMatrix<-", function(object, value) {
  standardGeneric("fullMatrix<-") 
})


#' qvalue object in experiment
#'
#' These generic functions access and set the \code{qvalue} object in the
#' \code{\linkS4class{edgeSet}} object.
#'
#' @usage qvalue.obj(object, ...)
#' 
#' @param object \code{\linkS4class{edgeSet}}
#' @param value S3 \code{object}: \code{\link{qvalue}} object
#' @param \dots Additional arguments used in \code{\link{qvalue}}.
#' 
#' @return  \code{qvalue.obj} returns a \code{\link{qvalue}} object.
#'     
#' @author John Storey, Andrew Bass
#' 
#' @seealso \code{\link{lrt}}, \code{\link{odp}} and 
#' \code{\linkS4class{edgeSet}}
#' 
#' @keywords qvalue.obj, qvalue.obj<-
#' 
#' @rdname qvalue.obj-method  
#' 
#' @exportMethod qvalue.obj                
setGeneric("qvalue.obj", function(object) standardGeneric("qvalue.obj"))

#' @rdname qvalue.obj-method
#' 
#' @exportMethod qvalue.obj<-
setGeneric("qvalue.obj<-", function(object, value) {
  standardGeneric("qvalue.obj<-") 
})                      

#' Individuals utilized in experiment
#'
#' These generic functions access and set the individual slot in 
#' \code{\linkS4class{edgeSet}}.
#'
#' @usage individual(object)
#' 
#' @param object \code{\linkS4class{edgeSet}}
#' @param value \code{factor}: individual identifier for each observation
#' 
#' @return \code{individual} returns information regarding individuals 
#' in the experiment.
#'    
#' @author John Storey, Andrew Bass
#' 
#' @seealso \code{\linkS4class{edgeSet}}
#' 
#' @keywords individual, individual<-
#' 
#' @rdname individual-method
#' 
#' @exportMethod individual
setGeneric("individual", function(object) standardGeneric("individual"))

#' @rdname individual-method
#' 
#' @exportMethod individual<-
setGeneric("individual<-", function(object, value) {
  standardGeneric("individual<-") 
})  
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