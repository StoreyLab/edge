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
  fitFull <- fitNull <- resNull <- resFull <- matrix(nrow=nr, ncol=n)
  for (i in 1:nr) {
    wlm_full <- lm.wfit(x = full.matrix, y = exprsData[i,], w = w[i,])
    wlm_null <- lm.wfit(x = null.matrix, y = exprsData[i,], w = w[i,])
    
    fitFull[i,] <- wlm_full$fitted.values
    fitNull[i,] <- wlm_null$fitted.values
    
    resFull[i,] <- wlm_full$residuals * sqrt(wlm_full$weights)
    resNull[i,] <- wlm_null$residuals * sqrt(wlm_full$weights)
  }
  dHFull <- diag(projMatrix(null.matrix))
  B.coef <- exprsData %*% full.matrix %*% ginv(t(full.matrix) %*% full.matrix)
  if (stat.var == "odp") {
    H.null <- projMatrix(null.matrix)
    full.matrix <- full.matrix - H.null %*% full.matrix
    full.matrix <- rm.zero.cols(full.matrix)
    H.full <- projMatrix(full.matrix)
    B.coef <- resNull %*% full.matrix %*% ginv(t(full.matrix) %*% full.matrix)
    dHFull <- diag(H.full)
    fitFull <- t(H.full %*% t(resNull))
    resFull <- resNull - fitFull
  }
  
  efObj <- new("deFit", fit.full = fitFull, fit.null = fitNull,
               dH.full = dHFull, res.full = resFull,
               res.null = resNull, beta.coef = B.coef,
               stat.type = stat.var)
  return(efObj)
}
