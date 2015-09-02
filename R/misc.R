bootstrap <- function(object, obs.fit, clustParms = NULL, bs.its = 100,
                      verbose = TRUE, mod.F = FALSE) {
  n.probes <- nrow(obs.fit@res.full)
  nf <- mod.df(object@full.matrix)
  null.stat <- matrix(nrow = n.probes,
                      ncol = bs.its)
  sType <- obs.fit@stat.type
  for (i in 1:bs.its) {
    if (verbose) {
      cat("\r", "Null iteration: ", i)
      if (i == bs.its) cat("\n")
    }
    exprs(object) <- null(obs.fit = obs.fit, nf = nf,
                          ind = object@individual)
    null.fit <- fit_models(object,
                           stat.type = sType)
    if (sType == "lrt") {
      post.var = NULL
      if (mod.F) {
        df_full <- ncol(object) - nf
        var_full <- rowSums(null.fit@res.full ^ 2) / df_full
        out <- squeezeVar(var_full, df_full, covariate = rowMeans(exprs(object)))
        post.var <- out$var.post
        prior.df <- out$df.prior
      }
      null.stat[, i] <- lrtStat(resNull = null.fit@res.null,
                                resFull = null.fit@res.full,
                                post.var = post.var)

    }
    else {
      null.stat[, i]  <- odpStat(n.res = null.fit@res.null,
                                 clustParms = clustParms)
    }
  }
  return(null.stat)
}
rescale <- function(x, sig) {
  means <- rowMeans(x)
  n <- ncol(x)
  rowsds <- sqrt((rowMeans(x ^ 2) - means ^ 2) * n / (n - 1))
  ret <- (x - means) * sig / rowsds + means
  return(ret)
}
null <- function(obs.fit, nf, ind) {
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
  df <- try(sum(diag(x %*% solve(t(x) %*% x) %*% t(x))), silent=TRUE)
  return(df)
}

createSet <- function(object, nMod=NULL, fMod=NULL, ind=NULL, grp=factor(NA)) {
  # Create deSet
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
  return(x[, colSums(abs(x)) > eps])
}


projMatrix <- function(x) {
  H <- x %*% ginv(t(x) %*% x) %*% t(x)
  H
}
