# Allows to set qvalue to S4 slot
setOldClass("qvalue")

deSetCheck <- function(object) {
  # Performs checks on an deSet object
  # 
  # Args:
  #   object: deSet object
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
  # Dimensionality test- this may be impossible to make in deSet
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

deFitCheck <- function(object) {
  # Performs checks on an deFit object
  # 
  # Args:
  #   object: deFit object
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

#' deSet class
#'
#' The deSet class was designed in order to complement the 
#' \code{\link{ExpressionSet}} class. While the \code{ExpressionSet} class 
#' contains information about the experiment, the deSet class 
#' contains both experimental information and additional information relevant 
#' to differential expression analysis.
#'
#' The deSet object is required for the following functions:
#' 
#'  @slot null.model \code{formula}: contains the adjustment variables in the 
#'  experiment.
#'  @slot full.model \code{formula}: contains the adjustment variables and the 
#'  biological variables of interest.
#'  @slot null.matrix \code{matrix}: the null model as a matrix.
#'  @slot full.matrix \code{matrix}: the full model as a matrix.
#'  @slot individual \code{factor}: containing information on individuals 
#'  sampled in the experiment.
#'  @slot qvalueObj S3 class \code{qvalue}: containing qvalue object. 
#'  See \code{\link{qvalue}} for additional details.
#'  
#' @section Methods:
#'  \describe{
#'  \item{\code{as(ExpressionSet, "deSet")}}{Coerce objects of ExpressionSet 
#'  to deSet}
#'  \item{\code{lrt(deSet, ...)}}{Likelihood ratio test}
#'  \item{\code{odp(deSet, ...)}}{Optimal discovery procedure}
#'  \item{\code{kl_clust(deSet, ...)}}{Clustering parameters for modular 
#'  optimal discovery procedure (mODP) method}
#'  \item{\code{fit_models(deSet, ...)}}{Linear regression for genes based on
#'              null and full models}
#'  \item{\code{apply_qvalue(deSet, ...)}}{Implements qvalue function on 
#'  deSet object}
#'  \item{\code{apply_snm(deSet, ...)}}{Implement surpervised normalization of
#'   microarrays on gene expression matrix on deSet object.}
#'  \item{\code{apply_sva(deSet, ...)}}{Estimate surrogate variables and adds 
#'  them to null/full models in a deSet object}
#'  \item{\code{fullMatrix(deSet)}}{Access and set full matrix from 
#'  deSet object}
#'  \item{\code{nullMatrix(deSet)}}{Access and set null matrix from 
#'  deSet object}
#'  \item{\code{fullModel(deSet)}}{Access and set full model from 
#'  deSet object}
#'  \item{\code{nullModel(deSet)}}{Access and set null model from 
#'  deSet object}
#'  \item{\code{individual(deSet)}}{Set individual slot from deSet object}
#'  \item{\code{qvalueObj(deSet)}}{Access qvalue object from deSet object. 
#'  See \code{\link{qvalue}}.}
#'  \item{\code{validObject(deSet)}}{Check validity of deSet object.}    
#'  }
#'  
#' @note 
#' The format for the model inputs are
#' \itemize{
#'  \item full.model: adjustment variables + biological variables 
#'  \item null.model: adjustment variables
#' }
#' The deSet object is created by either using \code{\link{deSet}}, 
#' \code{\link{build_models}}, or \code{\link{build_study}} functions. 
#' The qvalueObj is the slot of interest and is determined by either 
#' using the \code{odp} or the \code{lrt} function.
#' 
#' @author
#' John Storey, Jeffrey Leek, Andrew Bass
#' 
#' @seealso 
#' \code{\link{deSet}}
#' @inheritParams ExpressionSet
#' @exportClass deSet
setClass("deSet", slots=c(null.model = "formula", 
                          full.model = "formula",
                          null.matrix = "matrix",
                          full.matrix = "matrix",
                          individual = "factor", 
                          qvalueObj = "qvalue"),
         prototype=prototype(null.model = formula(NULL),
                             full.model = formula(NULL),
                             null.matrix = matrix(),
                             full.matrix = matrix(),
                             individual = as.factor(NULL),
                             qvalueObj = structure(list(), 
                                                    class = "qvalue")),
         validity = deSetCheck,
         contains = c("ExpressionSet"))

#' deFit class
#'
#' Object returned from \code{\link{fit_models}} function containing information 
#' regarding the model fits for the experiment. A least-squares algorithm 
#' is fit to both the full and null models of the experiment. 
#'
#' @section Slots: 
#'  \describe{
#'    \item{\code{fit.full}:}{Matrix containing fitted values for full model.}
#'    \item{\code{fit.null}:}{Matrix containing fitted values for null model.}
#'    \item{\code{res.full}:}{Matrix containing residuals for full model.}
#'    \item{\code{res.null}:}{Matrix containing residuals for null model.}
#'    \item{\code{dH.full}:}{Vector containing diagonal elements in projection 
#'    matrix for the full model.}
#'    \item{\code{beta.coef}:}{Matrix containing linear fitted coefficients 
#'    for full model.}
#'    \item{\code{stat.type}:}{String containing information on the statistic 
#'    of interest. Currently, the only options are ``lrt'' and ``odp''.}
#'  }
#'  
#' @section Methods:
#'  \describe{
#'  \item{\code{fitNull(deFit)}}{Access fitted data from null model}
#'  \item{\code{fitFull(deFit)}}{Access fitted data from full model}
#'  \item{\code{resNull(deFit)}}{Access null residuals from null model}
#'  \item{\code{resFull(deFit)}}{Access full residuals from full model}
#'  \item{\code{betaCoef(deFit)}}{Access beta coefficients in linear model}
#'  \item{\code{sType(deFit)}}{Access statistic type for model fitting utilized 
#'  in function}
#'  }
#' 
#' @note 
#' The deFit object can be used as additional inputs for the \code{lrt}, 
#' \code{odp} and \code{kl_clust} function but the object is automatically 
#' generated if not specified. 
#' 
#' @author 
#' John Storey, Jeffrey Leek, Andrew Bass
#' 
#' @seealso 
#' \code{\link{build_models}}
#' 
#' @exportClass deFit
setClass("deFit", slots=c(fit.full = "matrix", 
                          fit.null = "matrix", 
                          res.full = "matrix", 
                          res.null = "matrix",
                          dH.full = "vector",
                          beta.coef = "matrix",
                          stat.type = "character"),
         validity = deFitCheck)
