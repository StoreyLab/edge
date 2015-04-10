#' Formulates the experimental models
#'
#' \code{build_study} generates the full and null models for users unfamiliar 
#' with their experimental design. There are two types of experimental designs: 
#' static and time-course. For more details, refer to the user manual.  
#' 
#' @param data \code{matrix}: gene expression data.
#' @param sampling \code{string}: type of experiment. Either "static" or 
#' "timecourse". Default is "static".
#' @param grp \code{vector}: groups or biological variable in experiment. Optional.
#' @param tme \code{vector}: covariate of interest in time course study. Optional. 
#' @param ind \code{factor}: individual factor for repeated observations of 
#' the same individuals. Optional. 
#' @param bio.var \code{matrix}: biological variables. Optional.
#' @param basis.df \code{numeric}: degree of freedom of the spline fit for time 
#' course study. Default is 2.
#' @param basis.type \code{string}: either "ncs" or "ps" basis for time course
#'  study. Default is "ncs".
#' @param adj.var \code{matrix}: adjustment variables. Optional.
#'  
#' @return \code{\linkS4class{deSet}} object 
#'  
#' @examples 
#' # create ExpressionSet object from kidney dataset 
#' library(splines) 
#' data(kidney)
#' age <- kidney$age
#' sex <- kidney$sex
#' kidexpr <- kidney$kidexpr
#' 
#' # create deSet object from data
#' de_obj <- build_study(data = kidexpr, adj.var = sex, tme = age, 
#' sampling = "timecourse", basis.df = 4)
#' @seealso \code{\linkS4class{deSet}}, \code{\link{build_models}}
#' @author John Storey, Andy Bass 
#' @export
build_study = function(data, grp = NULL, adj.var = NULL, bio.var = NULL, 
                     tme = NULL, ind = NULL, 
                     sampling = c("static", "timecourse"), basis.df = 2, 
                     basis.type = c("ncs", "poly")) {
  n <- ncol(data)
  m <- nrow(data)
  if (!is.matrix(data)) {
    stop("data must be a matrix")
  }
  if (!is.null(tme)) {
    if (is.matrix(tme) | is.vector(tme)) {
      tme <- data.frame(tme)
    } else {
      stop("tme must be a matrix")
    }
   # intercept <- !apply(tme, 2, var)
   # tme <- subset(tme, select=!intercept)
  }
  if (!is.null(adj.var)) {
    if (is.matrix(adj.var) | is.vector(adj.var) | is.factor(adj.var)) {
      adj.var <- data.frame(adj.var)
    } else {
      stop("adj.var must be a matrix")
    }
    #intercept <- !apply(adj.var, 2, var)
   # adj.var <- subset(adj.var, select=!intercept)
  }
  if (!is.null(bio.var)) {
#    sampling <- "notApplicable"
    if (is.matrix(bio.var)| is.vector(bio.var) | is.factor(bio.var)) {
      bio.var <- data.frame(bio.var)
    } else {
      stop("bio.var must be a matrix")
    }
    #intercept <- !apply(bio.var, 2, var)
   # bio.var <- subset(bio.var, select=!intercept)
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
        grp <- data.frame(grp = as.factor(grp))
      } else {
        stop("grp must be a factor")
      }
     # intercept <- !apply(grp, 2, var)
      #grp <- subset(grp, select=!intercept)
    } else {
      if(sampling == "static") {
        stop("grp variable cannot be missing for static sampling.")
      }
      grp <- data.frame(grp=rep(1,n))
    }
    g <- nrow(unique(grp))
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
      if (length(varName) != 1) stop("Only one time variable is allowed. See ?deSet for information on how to create complicated models")
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
  exp_set <- new("ExpressionSet")
  pData(exp_set) <- data.frame(pdat)
  exprs(exp_set) <- as.matrix(data)
  edgeObj <- deSet(exp_set, full.model=as.formula(fmod), 
                     null.model=as.formula(nmod), individual=ind)
  return(edgeObj)  
}

#' Generate a deSet object with full and null models
#'
#' \code{build_models} creates a \code{\link{deSet}} object. The user inputs 
#' the full and null models. 
#' 
#' @param data \code{matrix}: gene expression data.
#' @param cov \code{data.frame}: the covariates in the study.
#' @param full.model \code{formula}: the adjustment and the biological 
#' variables of interest.
#' @param null.model \code{formula}: the adjustment variables. 
#' @param ind \code{factor}: individuals sampled in the study. Default is 
#' NULL. Optional.
#' 
#' @return \code{\linkS4class{deSet}} object 
#'  
#' @examples 
#' # create ExpressionSet object from kidney dataset
#' library(splines)
#' data(kidney)
#' age <- kidney$age
#' sex <- kidney$sex
#' kidexpr <- kidney$kidexpr
#' cov <- data.frame(sex = sex, age = age) 
#' 
#' # create models
#' null.model <- ~sex 
#' full.model <- ~sex + ns(age, df=4)
#' 
#' # create deSet object from data
#' de_obj <- build_models(data = kidexpr, cov = cov, null.model = null.model, 
#' full.model = full.model)
#' @seealso \code{\linkS4class{deSet}}, \code{\link{build_study}}
#' @author John Storey, Andy Bass 
#' @export
build_models <- function(data, cov, full.model = NULL, null.model = NULL, 
                      ind = NULL) {
  n <- ncol(data)
  m <- nrow(data)
  if (!is.matrix(data)) {
    stop("data must be a matrix")
  } else if (!is.data.frame(cov)) {
    stop("cov must be a data frame")
  } else if (is.null(full.model)) {
    stop("need an alternative model")
  }
  if (is.null(null.model)) {
    null.model <- ~1
  }
  if (!is(full.model, "formula") | !is(null.model, "formula")) {
    stop("alternative and null models must be formatted as a formula")
  }

  exp_set <- new("ExpressionSet")
  exprs(exp_set) <- data
  pData(exp_set) <- cov

  edgeObj <- deSet(exp_set, full.model = full.model, null.model = null.model, 
                   individual = ind)
  return(edgeObj)  
}
