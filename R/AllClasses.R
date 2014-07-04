# Allows to set qvalue to S4 slot
setOldClass("qvalue")

edgeSetCheck <- function(object) {
  # Performs checks on an edgeSet object
  # 
  # Args:
  #   object: edgeSet object
  #
  # Returns:
  #   TRUE/FALSE
  errors <- character()
  epsilon <- 10e-8
  # Allow easy conversion for an ExpressionSet using function 'as'
  if (is.list(object@null.model) && is.list(object@full.model) && 
        length(object@individual) == 0) {
    return(TRUE)
  }
  # Name mismatch
  f.vars <- all.vars(object@full.model)
  n.vars <- all.vars(object@null.model)
  names <- unique(c(f.vars, n.vars))
  if (sum((f.vars %in% c("grp", "bio.var", "time.basis"))) == 0) {
    if (sum(!(names %in% varLabels(object))) != 0) {
      msg <- paste("naming mismatch between phenoData covariates and models.")
      errors <- c(errors, msg) 
      return(errors)
    }
  }
  # Singular matrix
  xx0 <- model.matrix(object@null.model, data=object)
  xx1 <- model.matrix(object@full.model, data=object)
  #  sCheck.null <- min(svd(xx0)$d) < epsilon
  sCheck.full <- min(svd(xx1)$d) < epsilon
  #  if (sCheck.null) {
  #   msg <- paste("null model matrix is near singular.")
  #   errors <- c(errors, msg)
  #  }
  if (sCheck.full) {
    msg <- paste("full model matrix is near singular.")
    errors <- c(errors, msg)
  }
  # Dimensionality test- this may be impossible to make in edgeSet
  dataDim <- dim(exprs(object))
  if (dataDim[2] != nrow(xx1)) {
    msg <- paste( "dimension mismatch between full model and assayData.")
    errors <- c(errors, msg)
  }
  if (dataDim[2] != nrow(xx0)) {
    msg <- paste( "dimension mismatch between null model and assayData.")
    errors <- c(errors, msg)
  } 
  # inidividual input size
  if (length(object@individual) != 0) {
    if (length(object@individual) != ncol(exprs(object))) {
      msg <- paste("individual must be the same length as the number of arrays")
      errors <- c(errors, msg)
    }
  }
  if (length(errors) == 0) {
    TRUE
  } else {
    errors
  }
}

edgeFitCheck <- function(object) {
  # Performs checks on an edgeFit object
  # 
  # Args:
  #   object: edgeFit object
  #
  # Returns:
  #   TRUE/FALSE
  errors <- character()
  # Dimensionality test
  if (!(    (ncol(object@fit.full)==ncol(object@fit.null)
             && (ncol(object@res.full) == ncol(object@res.null))
             && (length(object@dH.full) == ncol(object@fit.full))
             && (ncol(object@fit.full) == ncol(object@res.null))))) {
    msg <- paste("column length of fitted matrices, dH.full and residuals",
                 "must be the same.")
    errors <- c(errors, msg) 
  }
  if (!((nrow(object@fit.full) == nrow(object@fit.null))
        && (nrow(object@res.full) == nrow(object@res.null))
        && (nrow(object@res.full) == nrow(object@fit.full))
        && (nrow(object@beta.coef) == nrow(object@fit.null)))) {
    msg <- paste("row length of fitted matrices and residuals", 
                 "must be the same.")
    errors <- c(errors, msg) 
  }
  # Correct statistic input check
  if (!(object@stat.type %in% c("lrt", "odp"))) {
    msg <- paste("stat.type must be lrt or odp. Inputted stat.type: ",
                 object@stat.type)
    errors <- c(errors, msg)
  }
  if (length(errors) == 0) {
    TRUE
  } else {
    errors
  }
}

#' edgeSet class
#'
#' Object inherits slots from an \code{ExpressionSet} object and has additional slots for 
#' the analysis steps. 
#'
#' The edgeSet object is the input for the main functions in the package edge. 
#' The format for the model inputs are
#' \itemize{
#'  \item full.model: adjustment variables + biological variables 
#'  \item null.model: adjustment variables
#' }
#' The \code{full.matrix} and \code{null.matrix} are created by the function
#' \code{\link{edgeSet}}, and the qvalue.obj is the slot of interest determined 
#' by either \code{odp} or \code{lrt} function.  
#'
#' @section Slots: 
#'  \describe{
#'    \item{\code{null.model}:}{\code{formula}: 
#'                              containing null model.}
#'    \item{\code{full.model}:}{\code{formula}:
#'                              containing full.model.}
#'    \item{\code{null.matrix}:}{\code{matrix}: 
#'                               containing null model data.}
#'    \item{\code{full.matrix}:}{\code{matrix}:
#'                               containing full model data.}
#'    \item{\code{individual}:}{\code{factor}: 
#'                              containing information on individuals in experiment.}
#'    \item{\code{qvalue.obj}:}{S3 class \code{qvalue}: 
#'                              containing qvalue object. See 
#'                              \code{\link{qvalue}}.}
#'    \item{\code{ExpressionSet}:}{Additional slots inherited from
#'                                 \code{\link{ExpressionSet}}.}                             
#'  }
#'  
#' @section Methods:
#'  \describe{
#'  \item{\code{as(exprSet, "edgeSet")}}{Coerce objects of ExpressionSet to
#'    edgeSet}
#'  \item{\code{lrt(edgeSet, ...)}}{Likelihood ratio test on edgeSet object}
#'  \item{\code{odp(edgeSet, ...)}}{Optimal discovery procedure on edgeSet object}
#'  \item{\code{klClust(edgeSet, ...)}}{Clustering parameters for modular optimal discovery procedure (mODP) method}
#'  \item{\code{edgeFit(edgeSet, ...)}}{Linear regression for genes based on
#'              null and full models}
#'  \item{\code{edgeQvalue(edgeSet, ...)}}{Implements qvalue function on edgeSet object}
#'  \item{\code{edgeSVA(edgeSet, ...)}}{Implement surpervised normalization of microarrays on gene expression
#'                                  matrix on edgeSet object.}
#'  \item{\code{edgeSNM(edgeSet, ...)}}{Estimate surrogate variables and adds them to null/full models
#'                                 from edgeSet object}
#'  \item{\code{fullMatrix(edgeSet)}}{Access and set full matrix from 
#'              edgeSet object}
#'  \item{\code{nullMatrix(edgeSet)}}{Access and set null matrix from 
#'              edgeSet object}
#'  \item{\code{fullModel(edgeSet)}}{Access and set full model from 
#'              edgeSet object}
#'  \item{\code{nullModel(edgeSet)}}{Access and set null model from
#'              edgeSet object}
#'  \item{\code{models(edgeSet)}}{Access null and full models from
#'              edgeSet object}
#'  \item{\code{individual(edgeSet)}}{Set individual slot from edgeSet object}
#'  \item{\code{dfModel(edgeSet)}}{Access the degree of freedom for null and full models from
#'              edgeSet object}
#'  \item{\code{qvalue.obj(edgeSet)}}{Access qvalue object from 
#'              edgeSet object. See \code{\link{qvalue}}.}
#'  \item{\code{validObject(edgeSet)}}{Check validity of edgeSet object.}    
#'  }
#'  
#' @note 
#' To create an edgeSet object, it is recommended to use the 
#' \code{\link{edgeSet}} function where the input is an ExpressionSet object, full 
#' and null models and optionally an individual variable. The edgeSet function generates
#' the full and null matrices. Once the edgeSet object is created, the user can
#' either use the \code{lrt} or \code{odp} method to determine the qvalue.obj 
#' slot in the edgeSet object. 
#' 
#' @author
#' John Storey, Jeffrey T. Leek, Andrew Bass
#' 
#' @seealso 
#' \code{\link{edgeSet}}
#' 
#' @keywords
#' edgeSet-class
#'
#' @rdname 
#' edgeSet-class
#' 
#' @exportClass
#' edgeSet 
setClass("edgeSet", slots=c(null.model = "formula", 
                            full.model = "formula",
                            null.matrix = "matrix",
                            full.matrix = "matrix",
                            individual = "factor", 
                            qvalue.obj = "qvalue"),
         prototype=prototype(null.model = formula(NULL),
                             full.model = formula(NULL),
                             null.matrix = matrix(),
                             full.matrix = matrix(),
                             individual = as.factor(NULL),
                             qvalue.obj = structure(list(), 
                                                    class = "qvalue")),
         validity = edgeSetCheck,
         contains = c("ExpressionSet"))

#' edgeFit class
#'
#' Object returned from \code{edgeFit-methods} function containing information 
#' on the model fits for the experiment. 
#'
#' Object contains output from fitting a linear model to each gene for both the 
#' full and null models. The stat.type caneither be "odp" or "lrt" depending on 
#' which statistic is of interest to generate the p-values.
#'
#' @section Slots: 
#'  \describe{
#'    \item{\code{fit.full}:}{Matrix of class \code{"matrix"}, 
#'                            containing fitted values for full model.}
#'    \item{\code{fit.null}:}{Matrix of class \code{"matrix"}, 
#'                            containing fitted values for null model.}
#'    \item{\code{res.full}:}{Matrix of class \code{"matrix"}, 
#'                            containing residuals for full model.}
#'    \item{\code{res.null}:}{Matrix of class \code{"matrix"}, 
#'                            containing residuals for null model.}
#'    \item{\code{dH.full}:}{Vector of class \code{"vector"}, 
#'                           containing diagonal elements in projection matrix 
#'                           for the full model.}
#'    \item{\code{beta.coef}:}{Matrix of class \code{"matrix"}, 
#'                           containing linear fitted coefficients for full model.}
#'    \item{\code{stat.type}:}{String of class \code{"string"}, 
#'                             containing information on the statistic of 
#'                             interest. Currently, the only options are lrt 
#'                             and odp.}
#'  }
#'  
#' @section Methods:
#'  \describe{
#'  \item{\code{modelFits(edgeFit, ...)}}{Fitted models of edgeFit object}
#'  \item{\code{fitNull(edgeFit)}}{Fitted data from null model}
#'  \item{\code{fitFull(edgeFit)}}{Fitted data from full model}
#'  \item{\code{resNull(edgeFit)}}{Null residuals from null model}
#'  \item{\code{resFull(edgeFit)}}{Full residuals from full model}
#'  \item{\code{betaCoef(edgeFit)}}{Beta coefficients in linear model}
#'  \item{\code{sType(edgeFit)}}{Statistic type for model fitting utilized in function}
#'  }
#' 
#' @note 
#' The edgeFit object can be the inputs for the \code{lrt}, \code{odp} 
#' and \code{klClust} function but the object is automatically generated if not
#' specified. 
#' 
#' @author 
#' John Storey, Jeffrey T. Leek, Andrew Bass
#' 
#' @seealso 
#' \code{\link{edgeFit}}
#' 
#' @keywords 
#' edgeFit-class
#' 
#' @rdname 
#' edgeFit-class 
#' 
#' @exportClass
#' edgeFit
setClass("edgeFit", slots=c(fit.full = "matrix", 
                            fit.null = "matrix", 
                            res.full = "matrix", 
                            res.null = "matrix",
                            dH.full = "vector",
                            beta.coef = "matrix",
                            stat.type = "character"),
         validity = edgeFitCheck)