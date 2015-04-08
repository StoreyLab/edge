# Allows to set qvalue to S4 slot
setOldClass("qvalue")

deSetCheck <- function(object) {
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

#' The differential expression class (deSet)
#'
#' The \code{deSet} class is designed in order to complement the 
#' \code{\link{ExpressionSet}} class. While the \code{ExpressionSet} class 
#' contains information about the experiment, the \code{deSet} class 
#' contains both experimental information and additional information relevant 
#' for differential expression analysis.
#' 
#' @slot null.model \code{formula}: contains the adjustment variables in the 
#' experiment. The null model is used for comparison when fitting the 
#' full model.
#' @slot full.model \code{formula}: contains the adjustment variables and the 
#' biological variables of interest.
#' @slot null.matrix \code{matrix}: the null model as a matrix.
#' @slot full.matrix \code{matrix}: the full model as a matrix.
#' @slot individual \code{factor}: contains information on which sample 
#' is from which individual in the experiment.
#' @slot qvalueObj \code{S3 object}: containing \code{qvalue} object. 
#' See \code{\link{qvalue}} for additional details.
#'  
#' @section Methods:
#'  \describe{
#'  \item{\code{as(ExpressionSet, "deSet")}}{Coerce objects of 
#'  \code{ExpressionSet} to \code{deSet}.}
#'  \item{\code{lrt(deSet, ...)}}{Performs a generalized likelihood ratio test 
#'  using the full and null models.}
#'  \item{\code{odp(deSet, ...)}}{Performs the optimal discovery procedure, 
#'  which is a new approach for optimally performing many hypothesis tests in 
#'  a high-dimensional study.}
#'  \item{\code{kl_clust(deSet, ...)}}{An implementation of mODP that assigns 
#'  genes to modules based off of the Kullback-Leibler distance.}
#'  \item{\code{fit_models(deSet, ...)}}{Fits a linear model to each gene by 
#'  method of least squares.}
#'  \item{\code{apply_qvalue(deSet, ...)}}{Applies \code{\link{qvalue}} 
#'  function.}
#'  \item{\code{apply_snm(deSet, ...)}}{Applies surpervised normalization of
#'   microarrays (\code{\link{snm}}) on gene expression data.}
#'  \item{\code{apply_sva(deSet, ...)}}{Applies surrogate variable analysis 
#'  (\code{\link{sva}}).}
#'  \item{\code{fullMatrix(deSet)}}{Access and set full matrix.}
#'  \item{\code{nullMatrix(deSet)}}{Access and set null matrix.}
#'  \item{\code{fullModel(deSet)}}{Access and set full model.}
#'  \item{\code{nullModel(deSet)}}{Access and set null model.}
#'  \item{\code{individual(deSet)}}{Access and set individual slot.}
#'  \item{\code{qvalueObj(deSet)}}{Access \code{qvalue} object. 
#'  See \code{\link{qvalue}}.}
#'  \item{\code{validObject(deSet)}}{Check validity of \code{deSet} object.}    
#'  }
#' 
#' @note
#' See \code{\link{ExpressionSet}} for slot information.
#' 
#' @author
#' John Storey, Jeffrey Leek, Andrew Bass
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

#' The differential expression class for the model fits
#'
#' Object returned from \code{\link{fit_models}} containing information 
#' regarding the model fits for the experiment.
#'
#' @slot fit.full \code{matrix}: containing fitted values for the full model.
#' @slot fit.null \code{matrix}: containing fitted values for the null model.
#' @slot res.full \code{matrix}: the residuals of the full model.
#' @slot res.null \code{matrix}: the residuals of the null model.
#' @slot dH.full \code{vector}: contains diagonal elements in the projection 
#' matrix for the full model.
#' @slot beta.coef \code{matrix}: fitted coefficients for the full model.
#' @slot stat.type \code{string}: information on the statistic of interest. 
#' Currently, the only options are ``lrt'' and ``odp''.
#'  
#' @section Methods:
#'  \describe{
#'  \item{\code{fitNull(deFit)}}{Access fitted data from null model.}
#'  \item{\code{fitFull(deFit)}}{Access fitted data from full model.}
#'  \item{\code{resNull(deFit)}}{Access residuals from null model fit.}
#'  \item{\code{resFull(deFit)}}{Access residuals from full model fit.}
#'  \item{\code{betaCoef(deFit)}}{Access beta coefficients in linear model.}
#'  \item{\code{sType(deFit)}}{Access statistic type of model fitting utilized 
#'  in function.}
#'  }
#' 
#' @author 
#' John Storey, Jeffrey Leek, Andrew Bass
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
