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
      #grp <- as.factor(grp)
      grp <- as.matrix(grp)
    } else {
      if(sampling == "static") {
        stop("grp variable cannot be missing for static sampling.")
      }
      grp <- as.matrix(rep(1, n), nrow=1)
    }
    g <- nrow(unique(grp))
    
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
#  intercept <- apply(as.matrix(adj.var), 2, function(x) length(unique(x)) == 1)
#  adj.var <- adj.var[, !intercept]
  if (sampling == "static") {
    if (g==1) {
      stop("grp must have more than one unique value for static sampling.")
    }
    # Create models for static
    pdat <- data.frame(adj.var, grp)
    if (is.null(adj.var)) {
      nmod <- ~1 
      fmod <- paste("~", paste(names(grp), collapse=" + "))
    } else {   
      fmod <- paste("~", paste(names(pdat), collapse=" + "))
      nmod <- paste("~", paste(names(adj.var), collapse=" + ")) 
    }
  }
  #need to get basis.df if it is NULL
  if (sampling == "timecourse") {
    basis.type <- match.arg(basis.type)
    if (basis.type == "ncs") {
      time.basis <- ns(tme, df=, intercept=FALSE)
      time.basis <- paste("ns(tme, df=", basis.df,", intercept=FALSE)", sep="")
    } else if (basis.type == "poly") {
      time.basis <- paste("bs(tme, df=", basis.df,", intercept=FALSE)", sep="")
    }
 #   tName <- NULL
  #  for (i in 1:basis.df){
  #    tName <- c(tName, paste("time.basis.", i, sep=""))
  #  }
  #  colnames(time.basis) <- tName
    
    if (g == 1) {
      # time course with no groups
      pdat <- data.frame(adj.var,  tme) 
      if (is.null(adj.var)) {
        nmod <- "~1" 
        fmod <- paste("~", time.basis) 
      } else {
        fmod <- paste("~", paste(names(adj.var), collapse=" + "), "+", time.basis) 
        nmod <- paste("~", paste(names(adj.var), collapse=" + ")) 
      }  
    } else {
      if (is.null(adj.var)) {
        pdat <- data.frame(grp)
      } else {
        pdat <- data.frame(adj.var, grp) 
      }
      # time course with groups
      nmod <- paste(paste("~", paste(names(pdat), collapse=" + ")), "+", time.basis)
      fmod <- paste(paste("~", paste(names(pdat), collapse=" + ")),"+",time.basis,"+", paste( "(", paste(names(data.frame(grp)), collapse=" + "), ")", ":", time.basis))  }
      pdat$tme <- tme
  }
  expSet <- new("ExpressionSet")
  pData(expSet) <- data.frame(pdat)
  exprs(expSet) <- as.matrix(data)
  edgeObj <- edgeSet(expSet, full.model=as.formula(fmod), null.model=as.formula(nmod), individual=(ind))
  return(edgeObj)  
}