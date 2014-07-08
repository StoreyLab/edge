#' Creates an edgeSet object and formulates appropriate models for users
#'
#' \code{edgeModel} is a function to create an edge object without an ExpressionSet. Alternative and null models are created based on experiment type: Either "static" or "timecourse". For more detail refer to the user manual.  
#' 
#' @param data matrix- Gene expression data.
#' @param sampling string- Type of experiment. Either "static" or "timecourse". Default is "static".
#' @param grp vector- Grouping variable in experiment. Optional.
#' @param tme vector- Covariate of interest in time course study. Optional. 
#' @param ind factor- Individual factor. Optional. 
#' @param basis.df numeric- Degree of freedom of the spline fit for time course study. Default is 2.
#' @param basis.type string- Either "ncs" or "ps" basis for time course study. Default is "ncs".
#' @param bio.var matrix- Biological variables. Optional.
#' @param adj.var matrix- Adjustment Variables. Optional.
#'  
#' @return \code{edgeModel} returns an \code{\linkS4class{edgeSet}} object with the following slots assigned:
#'   \describe{
#'     \item{\code{full.model:}}{alternative model equation}
#'    \item{\code{null.model:}}{null model equation}
#'    \item{\code{full.matrix:}}{alternative model in matrix form}
#'    \item{\code{null.matirx:}}{null model in matrix form}
#'    \item{\code{individual:}}{individuals in experiment (factor)}
#'    \item{\code{ExpressionSet:}}{inherits ExpressionSet object (assayData, phenoData) created in function}
#'  }
#'  
#' @examples 
#' # Create ExpressionSet object from kidney dataset 
#' library(splines) 
#' data(kidney)
#' sex <- kidney$sex[kidney$tissue == "c"]
#' age <- kidney$age[kidney$tissue == "c"]
#' kidexpr <- log(kidney$kidexpr[, kidney$tissue == "c"] + 10)
#' 
#' #Create edgeSet object from data
#' edgeObj <- edgeModel(data=kidexpr, adj.var=model.matrix(~sex), tme=age, sampling="timecourse", basis.df=4)
#' @name edgeModel
#' @rdname edgeModel
#' @seealso \code{\link{edgeSet}}
#' @author John Storey, Andy Bass 
#' @aliases edgeModel
#' @export
edgeModel = function(data, sampling=c("static", "timecourse"), grp = NULL, tme=NULL, ind=NULL, basis.df=2, basis.type = c("ncs", "poly"), bio.var=NULL, adj.var=NULL)
{
  n <- ncol(data)
  m <- nrow(data)
  if (!is.null(bio.var)) {
    sampling <- "notApplicable"
    if (is.vector(bio.var)) {
      bio.var <- matrix(bio.var, ncol=1)
    }
    intercept <- apply(bio.var, 2, function(x) length(unique(x)) == 1)
    if (sum(intercept) > 0) {
      warning("Removing intercept column from bio.var")
    }
    bio.var <- bio.var[, !intercept]
    if (is.null(adj.var)) {
      adj.var <- matrix(rep(1, n), ncol=1)
    }
    if (is.vector(adj.var)) {
      adj.var = matrix(adj.var, ncol=1)
    }
    # Create models
    fmod <- ~adj.var + bio.var
    nmod <- ~adj.var
    pdat <- cbind(adj.var, bio.var)
  } else {
    sampling <- match.arg(sampling, choices=c("static", "timecourse"))
    if (!is.null(grp)) {
      grp <- as.factor(grp)
    } else {
      if(sampling == "static") {
        stop("grp variable cannot be missing for static sampling.")
      }
      grp <- rep(1, n)
    }
    g <- length(unique(grp))
    
  #  if (!is.null(adj.var)) {
  #    if (is.vector(adj.var)) {
  #      adj.var <- matrix(adj.var, ncol=1)
  #    }
  #    intercept <- as.logical(max(apply(adj.var, 2, function(x) length(unique(x))) == 1))
  #    if (!intercept) {
  #      adj.var <- cbind(rep(1, n), adj.var)
  #    }
  #  }
   # if (is.null(adj.var)) {
    #  adj.var <- matrix(rep(1, n), ncol=1)
   # }   
    #if (!is.null(ind)) {
    #  if (is.vector(ind)) {
    #    ind <- matrix(ind, ncol=1)
    #  }
    #  ind2 <- NULL
    #  for (i in 1:ncol(ind)) {
    #    ind2 <- cbind(ind2, model.matrix(~ -1 + as.factor(ind[, i])))
    #  }
     # ind <- ind2
    #} 
  }
  if (sampling == "static") {
    if (g==1) {
      stop("grp must have more than one unique value for static sampling.")
    }
    # Create models for static
    if (is.null(adj.var)) {
      nmod <- ~1 
      fmod <- ~grp
    } else {   
      fmod <- ~adj.var + grp
      nmod <- ~adj.var  
    }
    pdat <- cbind(adj.var, grp)
  }
  #need to get basis.df if it is NULL
  if (sampling == "timecourse") {
    basis.type <- match.arg(basis.type)
    if (basis.type == "ncs") {
      knts = quantile(tme, probs=seq(0, 1, length=(basis.df + 1)))[-c(1, (basis.df + 1))]
      time.basis <- ns(tme, knots=knts, intercept=FALSE)
    } else if (basis.type == "poly") {
      knts <- quantile(tme, probs=seq(0, 1, length=(basis.df + 1)))[-c(1, (basis.df + 1))]
      time.basis <- bs(tme, knots=knts, intercept=FALSE)
    }
    if (g == 1) {
      # time course with no groups
      if (is.null(adj.var)) {
        nmod <- ~1 
        fmod <- ~time.basis
      } else {
        fmod <- ~adj.var + time.basis
        nmod <- ~adj.var
      }  
      pdat <- cbind(adj.var,  time.basis) 
    } else {
      # time course with groups
      if (is.null(adj.var)) {
        nmod <- ~grp + time.basis 
        fmod <- ~grp + time.basis + time.basis:grp
      } else {
        nmod <- ~adj.var + grp + time.basis 
        fmod <- ~adj.var + grp + time.basis + time.basis:grp
      }
      pdat <- cbind(adj.var, as.matrix(grp), time.basis)   
    }
  }
  expSet <- new("ExpressionSet")
  pData(expSet) <- data.frame(pdat)
  exprs(expSet) <- as.matrix(data)
  edgeObj <- edgeSet(expSet, full.model=fmod, null.model=nmod, individual=(ind))
  return(edgeObj)  
}