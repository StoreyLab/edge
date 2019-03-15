fit_wmodels <- function(object, w = NULL,  stat.type = c("lrt", "odp")) {
  exprsData <- exprs(object)
  n <- ncol(exprsData)
  nr <- nrow(exprsData)
  stat.var <- match.arg(stat.type, c("lrt", "odp"))
  null.matrix <- object@null.matrix
  full.matrix <- object@full.matrix
  if (length(object@individual) != 0) {
    ind.matrix <- model.matrix(~-1 + as.factor(object@individual))
    Hi <- projMatrix(ind.matrix)
    fitInd <- t(Hi %*% t(exprsData))
    exprsData <- exprsData - fitInd
    full.matrix <- full.matrix - Hi %*% full.matrix
    null.matrix <- null.matrix - Hi %*% null.matrix
    full.matrix <- rm.zero.cols(full.matrix)
    null.matrix <- rm.zero.cols(null.matrix)
  }
  fitFull <- fitNull <- resNull <- resFull <- dHFull <- matrix(nrow=nr, ncol=n)
  for (i in 1:nr) {
    wlm_null <- lm.wfit(x = null.matrix, y = exprsData[i,], w = w[i,])
    fitNull[i,] <- wlm_null$fitted.values
    resNull[i,] <- wlm_null$residuals * sqrt(wlm_null$weights)
    if (stat.var != "odp") {
      wlm_full <- lm.wfit(x = full.matrix, y = exprsData[i,], w = w[i,])
      dHFull[i,] <- diag(projMatrix(sqrt(w[i,]) * full.matrix))# double check
      fitFull[i,] <- wlm_full$fitted.values
      B.coef <- matrix(NA, ncol = length(w[i,]))#wlm_full$coefficients
      resFull[i,] <- wlm_full$residuals * sqrt(wlm_full$weights)
    } else {
      # W <- diag(sqrt(w[i,]))
      w_sqrt <- sqrt(w[i,])
      f.matrix.scaled <- full.matrix * w_sqrt
      H.null <- projMatrix(null.matrix * w_sqrt)
      f.matrix.scaled <- f.matrix.scaled - H.null %*% f.matrix.scaled
      f.matrix.scaled <- rm.zero.cols(f.matrix.scaled)
      H.full <- projMatrix(f.matrix.scaled)
      res.n <- wlm_null$residuals * w_sqrt
      B.coef <- matrix(NA, ncol = length(w_sqrt))#res.n %*% full.matrix.scaled %*% ginv(t(full.matrix.scaled) %*% full.matrix.scaled)
      dHFull[i,] <- diag(H.full)
      fitFull[i,] <- H.full %*% res.n
      resFull[i,] <- res.n - fitFull[i,]
    }
  }
  efObj <- new("deFit", fit.full = fitFull, fit.null = fitNull,
               dH.full = dHFull, res.full = resFull,
               res.null = resNull, beta.coef = B.coef,
               stat.type = stat.var)
  return(efObj)
}
