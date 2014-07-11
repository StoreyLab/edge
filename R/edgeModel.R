#' Creates an edgeSet object and formulates appropriate models for users
#'
#' \code{edgeModel} is a function to create an edge object without an ExpressionSet. Alternative and null models are created based on experiment type: Either "static" or "timecourse". For more detail refer to the user manual.  
#' 
#' @param data matrix- Gene expression data.
#' @param sampling string- Type of experiment. Either "static" or "timecourse". Default is "static".
#' @param grp vector- Groups or biological variable in experiment. Optional.
#' @param tme vector- Covariate of interest in time course study. Optional. 
#' @param ind factor- Individual factor. Optional. 
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
#' edgeObj <- edgeModel(data=kidexpr, adj.var=model.matrix(~sex), tme=age, sampling="timecourse", basis.df=4)
#' @name edgeModel
#' @rdname edgeModel
#' @seealso \code{\link{edgeSet}}
#' @author John Storey, Andy Bass 
#' @aliases edgeModel
#' @export
edgeModel = function(data, sampling=c("static", "timecourse"), tme=NULL, ind=NULL, basis.df=2, basis.type = c("ncs", "poly"), grp=NULL, adj.var=NULL) {
  n <- ncol(data)
  m <- nrow(data)
  sampling <- match.arg(sampling, choices=c("static", "timecourse"))
  if (!is.null(grp)) {
    grp <- as.matrix(grp)
  } else {
    if(sampling == "static") {
      stop("grp variable cannot be missing for static sampling.")
    }
    grp <- as.matrix(rep(1, n), nrow=1)
  }
  g <- nrow(unique(grp))

  if (sampling == "static") {
    if (g==1) {
      stop("grp must have more than one unique value for static sampling.")
    }
    
   # pdat <- data.frame(adj.var, grp)
    if (is.null(adj.var)) {
      pdat <- data.frame(grp)
      nmod <- ~1 
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
      pdat <- data.frame(adj.var, tme) 
      if (is.null(adj.var)) {
        nmod <- "~1" 
        fmod <- paste("~", time.basis) 
      } else {
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
  expSet <- new("ExpressionSet")
  pData(expSet) <- data.frame(pdat)
  exprs(expSet) <- as.matrix(data)
  edgeObj <- edgeSet(expSet, full.model=as.formula(fmod), null.model=as.formula(nmod), individual=(ind))
  return(edgeObj)  
}