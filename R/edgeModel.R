#' Creates an edgeSet object and formulates appropriate models for users
#'
#' \code{edgeModel} is a function to create an edge object without an ExpressionSet. Alternative and null models are created based on experiment type: Either "static" or "timecourse". For more detail refer to the user manual.  
#' 
#' @param data matrix- Gene expression data.
#' @param sampling string- Type of experiment. Either "static" or "timecourse". Default is "static".
#' @param grp vector- Groups or biological variable in experiment. Optional.
#' @param tme vector- Covariate of interest in time course study. Optional. 
#' @param ind factor- Individual factor. Optional. 
#' @param bio.var matrix- Biological variables. Optional.
#' @param basis.df numeric- Degree of freedom of the spline fit for time course study. Default is 2.
#' @param basis.type string- Either "ncs" or "ps" basis for time course study. Default is "ncs".
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
#' sex <- kidney$sex
#' age <- kidney$age
#' kidexpr <- kidney$kidexpr
#' 
#' #Create edgeSet object from data
#' edgeObj <- edgeModel(data=kidexpr, adj.var=data.frame(sex), tme=age, sampling="timecourse", basis.df=4)
#' @name edgeModel
#' @rdname edgeModel
#' @seealso \code{\link{edgeSet}}
#' @author John Storey, Andy Bass 
#' @aliases edgeModel
#' @export
edgeModel = function(data, sampling=c("static", "timecourse"), tme=NULL, ind=NULL, basis.df=2, basis.type = c("ncs", "poly"), grp=NULL, adj.var=NULL, bio.var=NULL) {
  n <- ncol(data)
  m <- nrow(data)
  if (!is.matrix(data)) {
    stop("data must be a matrix")
  }
  if (!is.null(tme)) {
    if (is.matrix(tme) | is.vector(tme)) {
      tme <- data.frame(tme)
    } else {
      stop("tme must be a matrix/vector")
    }
    intercept <- !apply(tme, 2, var)
    tme <- subset(tme, select=!intercept)
  }
  
  if (!is.null(adj.var)) {
    if (is.matrix(adj.var) | is.data.frame(adj.var) | is.vector(adj.var)) {
      adj.var <- data.frame(adj.var)
    } else {
      stop("adj.var must be a matrix/vector/factor/data.frame")
    }
    intercept <- !apply(adj.var, 2, var)
    adj.var <- subset(adj.var, select=!intercept)
  }
  if (!is.null(bio.var)) {
    sampling <- "notApplicable"
    if (is.matrix(bio.var) | is.vector(bio.var)) {
      bio.var <- data.frame(bio.var)
    } else {
      stop("bio.var must be a matrix/vector/factor")
    }
    intercept <- !apply(bio.var, 2, var)
    bio.var <- subset(bio.var, select=!intercept)
    # Create models
    if (is.null(adj.var)) {
      pdat <- data.frame(bio.var)
      fmod <- paste("~", paste(names(pdat), collapse=" + "))
      nmod <- "~1"
    } else {
      pdat <- data.frame(adj.var, bio.var)
      fmod <- paste("~", paste(names(pdat), collapse=" + "))
      nmod <- paste("~", paste(names(adj.var), collapse=" + ")) 
    }
  } else {
    sampling <- match.arg(sampling, choices=c("static", "timecourse")) 
    if (!is.null(grp)) {
      if (is.factor(grp)) {
        grp <- data.frame(model.matrix(~grp))
      } else {
        stop("grp must be a factor")
      }
      intercept <- !apply(grp, 2, var)
      grp <- subset(grp, select=!intercept)
    } else {
      if(sampling == "static") {
        stop("grp variable cannot be missing for static sampling.")
      }
      grp <- data.frame(grp=rep(1,n))
    }
    g <- length(unique(grp))
    if (sampling == "static") {
      if (g==1) {
        stop("grp must have more than one unique value for static sampling.")
      }
      if (is.null(adj.var)) {
        pdat <- data.frame(grp)
        nmod <- "~1" 
        fmod <- paste("~", paste(names(pdat), collapse=" + "))
      } else {  
        pdat <- data.frame(adj.var, grp)
        fmod <- paste("~", paste(names(pdat), collapse=" + "))
        nmod <- paste("~", paste(names(adj.var), collapse=" + ")) 
      }
    }
  
    if (sampling == "timecourse") {
      basis.type <- match.arg(basis.type)
      varName <- colnames(data.frame(tme))
      if (length(varName) != 1) stop("Only one time variable is allowed. See ?edgeSet for information on how to create complicated models")
      if (basis.type == "ncs") {
        time.basis <- paste("ns(", varName,", df=", basis.df,", intercept=FALSE)", sep="")
      } else if (basis.type == "poly") {
        time.basis <- paste("bs(", varName,", df=", basis.df,", intercept=FALSE)", sep="")
      }
      if (g == 1) {
        # time course with no groups
        if (is.null(adj.var)) {
          pdat <- data.frame(tme) 
          nmod <- "~1" 
          fmod <- paste("~", time.basis) 
        } else {
          pdat <- data.frame(adj.var, tme) 
          fmod <- paste("~", paste(names(adj.var), collapse=" + "), "+", time.basis) 
          nmod <- paste("~", paste(names(adj.var), collapse=" + ")) 
        }  
      } else {
        if (is.null(adj.var)) {
          pdat <- data.frame(tme, grp)
        } else {
          pdat <- data.frame(tme, adj.var, grp) 
        }
        # time course with groups
        nmod <- paste(paste("~", paste(names(pdat)[-1], collapse=" + ")), "+", time.basis)
        fmod <- paste(paste("~", paste(names(pdat)[-1], collapse=" + ")),"+",time.basis,"+", paste( "(", paste(names(pdat)[ncol(pdat)], collapse=" + ", sep=""), ")", ":", time.basis))  }
    }
  }
  expSet <- new("ExpressionSet")
  pData(expSet) <- data.frame(pdat)
  exprs(expSet) <- as.matrix(data)
  edgeObj <- edgeSet(expSet, full.model=as.formula(fmod), null.model=as.formula(nmod), individual=ind)
  return(edgeObj)  
}