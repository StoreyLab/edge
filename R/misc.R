bootstrap <- function(object, obs.fit, clustParms = NULL, bs.its = 100, 
                      verbose = TRUE) {
  # lrt Null statistic
  #
  # Args:
  #   object: class edgeSet
  #   obs.fit: class edgeFit
  #   df: list- d.o.f of null and full models and residuals
  #   bs.its: numeric- number of iterations, default is 100
  #   verbose: string- if TRUE, print iteration number 
  #
  # Returns:
  #   matrix- null statistic 
  n.probes <- nrow(obs.fit@res.full)
  null.stat <- matrix(nrow = n.probes,
                      ncol = bs.its)
  nf <-  mod.df(object@full.matrix)
  sType <- obs.fit@stat.type
  for (i in 1:bs.its) {
    if (verbose) {
      cat("\r", "Null iteration: ", i)
      if (i == bs.its) cat("\n")
    }
    exprs(object) <- null(obs.fit = obs.fit, nf = nf,
                          ind = object@individual)
    null.fit <- edgeFit(object, 
                        stat.type = sType)
    if (sType == "lrt") {
      null.stat[, i]  <- lrtStat(resNull = null.fit@res.null, 
                                 resFull = null.fit@res.full)
    }
    else {
      null.stat[, i]  <- odpStat(n.res = null.fit@res.null, 
                                 clustParms = clustParms)
    }
  }
  return(null.stat)
}    
rescale <- function(x, sig) {
  # rescale residuals
  # 
  # Args:
  #   x: matrix- residuals
  #   sig: vector- sd
  #
  # Returns:
  #   matrix- rescaled residuals for ODP
  
  means <- rowMeans(x)
  n <- ncol(x)
  rowsds <- sqrt((rowMeans(x ^ 2) - means ^ 2) * n / (n - 1))
  ret <- (x - means) * sig / rowsds + means
  return(ret)
}
null <- function(obs.fit, nf, ind) {
  # Calculates null data
  #
  # Args:
  #   obs.fit: edgeFit class
  #   df.mod: list- degree of freedom of bio/adj matrices
  #   ind: factor- individuals 
  # 
  # Returns:
  #   null.dat: matrix- simulated null data 
  # Initializations  
  stat.var <- obs.fit@stat.type
  n <- ncol(obs.fit@res.full)
  if (sum(!is.na(ind[1])) > 0) {
    ind <- model.matrix(~-1 + as.factor(ind))
    wts <- sqrt(1 - diag(ind %*% solve(t(ind) %*% ind) %*% t(ind)))
  } else {
    ind <- NULL
    wts <- rep(1, n)
  }
  wts <- wts * sqrt(1 - obs.fit@dH.full) 
  res.full <- t(t(obs.fit@res.full) / wts)
  # Random mix columns of residuals from full model
  vv <- sample(1:n, replace = TRUE)
  bs.res <- res.full[, vv]
  # Add random residuals to null data
  if (stat.var == "lrt") {
    null.dat <- obs.fit@fit.null + bs.res
  } else {
    sig1 <- sqrt(rowSums(obs.fit@res.full ^ 2) / (n - nf))
    bs.res <- rescale(x = bs.res, 
                      sig = sig1)
    null.dat <- obs.fit@fit.null + bs.res
  }
  return(null.dat)
}

mod.df <- function(x) {
  # Degree of freedom
  #
  # Args:
  #   x: matrix
  #
  # Returns:
  #   df: numeric- degrees of freedom of x
  df <- try(sum(diag(x %*% solve(t(x) %*% x) %*% t(x))), silent=TRUE)
  return(df)
}

createSet <- function(object, nMod=NULL, fMod=NULL, ind=NULL, grp=factor(NA)) {
  # Create edgeSet
  #  require(splines)
  object@null.model <- nMod
  object@full.model <- fMod
  mmf <- model.matrix(object = fMod, data = object)
  mmn <- model.matrix(object = nMod, data = object)
  colnames(mmf) <- NULL
  colnames(mmn) <- NULL
  object@null.matrix <- mmn
  object@full.matrix <- mmf 
  object@individual <- as.factor(ind)
  validObject(object)
  object
}

rm.zero.cols <- function(x, eps = 10e-12) {
  # Remove null columns
  #
  # Args:
  #   x: matrix
  #   eps: numeric- threshold
  # 
  # Results:
  #   matrix- with zero columns removed
  return(x[, colSums(abs(x)) > eps])
}


projMatrix <- function(x) {
  H <- x %*% ginv(t(x) %*% x) %*% t(x)
  H
}