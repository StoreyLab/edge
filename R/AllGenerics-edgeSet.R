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