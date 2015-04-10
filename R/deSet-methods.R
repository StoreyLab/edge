#' @rdname fit_models
setMethod("fit_models",
          "deSet",
          function(object, stat.type = c("lrt", "odp")) {
            # Initializations
            stat.var <- match.arg(stat.type, c("lrt", "odp"))
            exprsData <- exprs(object)
            n <- ncol(exprsData)
            null.matrix <- object@null.matrix
            full.matrix <- object@full.matrix
            # Rescale if there are group individual factors
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
            # Fitted exprsData and statistics under null model and full model
            H.null <- projMatrix(null.matrix)
            fitNull <- t(H.null %*% t(exprsData))
            resNull <- exprsData - fitNull
            if (stat.var == "odp") {
              full.matrix <- full.matrix - H.null %*% full.matrix
              full.matrix <- rm.zero.cols(full.matrix)
              H.full <- projMatrix(full.matrix)
              B.coef <- resNull %*% full.matrix %*% ginv(t(full.matrix) %*% full.matrix)
              dHFull <- diag(H.full)
              fitFull <- t(H.full %*% t(resNull))
              resFull <- resNull - fitFull
            } else {
              H.full <- projMatrix(full.matrix)
              dHFull <- diag(H.full)
              B.coef <- exprsData %*% full.matrix %*% ginv(t(full.matrix) %*% full.matrix)
              fitFull <- t(H.full %*% t(exprsData))
              resFull <- exprsData - fitFull
            }
            efObj <- new("deFit", fit.full = fitFull, fit.null = fitNull,
                         dH.full = dHFull, res.full = resFull,
                         res.null = resNull, beta.coef = B.coef,
                         stat.type = stat.var)
            return(efObj)
          })

#' @rdname odp
setMethod("odp",
          signature = signature(object = "deSet", de.fit = "missing"),
          function(object, de.fit, odp.parms = NULL, bs.its = 100,
                   n.mods = 50, seed = NULL, verbose = TRUE, ...)  {
            de.fit <- fit_models(object,
                                 stat.type = "odp", ...)
            results <- odp(object, de.fit,
                           odp.parms = odp.parms,
                           n.mods = n.mods,
                           bs.its = bs.its,
                           seed = seed,
                           verbose = verbose, ...)
            return(results)
          })

#' @rdname odp
setMethod("odp",
          signature = signature(object = "deSet", de.fit = "deFit"),
          function(object, de.fit, odp.parms = NULL, bs.its = 100,
                   n.mods = 50, seed = NULL, verbose = TRUE, ...) {
            if (!is.null(seed)) {
              set.seed(seed)
            }
            if (is.null(odp.parms)) {
              odp.parms <- kl_clust(object,
                                   n.mods = n.mods)
            } else if (sum(!(names(odp.parms) %in% c("mu.full", "sig.full",
                                                     "mu.null", "sig.null",
                                                     "n.per.mod",
                                                     "clustMembers"))) != 0) {
              stop("Not a correct ODP parameter list. See kl_clust documentation")
            }
            odp.stat <- odpStat(n.res = de.fit@res.null,
                                clustParms = odp.parms)
            null.stat <- bootstrap(object = object,
                                   obs.fit = de.fit,
                                   clustParms = odp.parms,
                                   bs.its = bs.its,
                                   verbose = verbose)
            pval <- empPvals(stat = odp.stat,
                             stat0 = null.stat, ...)
            qvalueObj(object) <- qvalue(p = pval, ...)
            return(object)
          })

#' @rdname lrt
setMethod("lrt",
          signature = signature(object = "deSet", de.fit = "missing"),
          function(object, de.fit, nullDistn = c("normal", "bootstrap"),
                   bs.its = 100, seed = NULL, ...) {
            de.fit <- fit_models(object,
                                  stat.type = "lrt")
            results <- lrt(object,
                           de.fit = de.fit,
                           nullDistn = nullDistn,
                           bs.its = bs.its,
                           seed = seed,
                           verbose = verbose, ...)
            return(results)
          })

#' @rdname lrt
setMethod("lrt",
          signature = signature(object = "deSet", de.fit = "deFit"),
          function(object, de.fit, nullDistn = c("normal", "bootstrap"),
                   bs.its = 100, seed = NULL, verbose = TRUE, ...) {
            # Initilizations
            nFull <- ncol(object@full.matrix)
            nNull <- ncol(object@null.matrix)
            n <- ncol(object)
            m <- nrow(object)
            if (!is.null(seed)) {
              set.seed(seed)
            }
            nullDistn <- match.arg(nullDistn, c("normal", "bootstrap"))
            # lrt observed stat
            stat <- lrtStat(resNull = de.fit@res.null,
                            resFull = de.fit@res.full)
            # If nullDistn is normal then return p-values from F-test else
            # return empirical p-values from qvalue package
            if (nullDistn == "normal") {
              df1 <- nFull - nNull
              df2 <- n - nFull
              stat <- stat * df2 / df1
              pval <- 1 - pf(stat,
                             df1 = df1,
                             df2 = df2)
              qvalueObj(object) <- qvalue(p = pval, ...)
              return(object)
            } else {
              null.stat <- bootstrap(object = object,
                                     obs.fit = de.fit,
                                     bs.its = bs.its,
                                     verbose = verbose)
              pval <- empPvals(stat = stat,
                               stat0 = null.stat, ...)
              qvalueObj(object) <- qvalue(pval, ...)
              return(object)
            }
          })

#' @rdname kl_clust
setMethod("kl_clust",
          signature = signature(object = "deSet", de.fit = "missing"),
          function(object, de.fit,  n.mods = 50)  {
            de.fit <- fit_models(object, stat.type = "odp")
            results <- kl_clust(object, de.fit,
                               n.mods = n.mods)
            return(results)
          })

#' @rdname kl_clust
setMethod("kl_clust",
          signature = signature(object = "deSet", de.fit = "deFit"),
          function(object, de.fit, n.mods = 50) {
            nf <- mod.df(object@full.matrix)
            nn <- mod.df(object@null.matrix)
            mod.member <- klmod(de.fit, nf = nf,
                                n.mods = n.mods)
            return(mod.parms(de.fit, nf = nf, nn=nn,
                             clMembers = mod.member))
          })

#' @rdname summary
setMethod("summary",
          signature = signature(object="deSet"),
          function(object, ...) {
            cat('\n')
            cat('ExpressionSet Summary', '\n', '\n')
            tmp <- as(object, "ExpressionSet")
            print(tmp)
            cat('\n')
            cat('de Analysis Summary', '\n', '\n')
            cat('Total number of arrays:', ncol(exprs(object)), '\n')
            cat('Total number of probes:', nrow(exprs(object)), '\n', '\n')
            cat('Biological variables:', '\n')
            cat('\tNull Model:')
            print(nullModel(object))
            cat('\n\tFull Model:')
            print(fullModel(object))
            cat('\n')
            if (length(object@individual) != 0) {
              cat('Individuals:', '\n')
              ind <- as.numeric(object@individual)
              print(individual(object))
              cat('\n')
            }
            cat('.......', '\n', '\n')
            if (!is.null(object@qvalueObj$pvalues)) {
              cuts <- c(0.0001, 0.001, 0.01, 0.025, 0.05, 0.10, 1)
              digits <- getOption("digits")
              cat("\nStatistical significance summary:\n")
              cat("pi0:", format(object@qvalueObj$pi0, digits = digits),
                  "\n", sep = "\t")
              cat("\n")
              cat("Cumulative number of significant calls:\n")
              cat("\n")
              counts <- sapply(cuts, function(x) c("p-value" = sum(object@qvalueObj$pvalues < x),
                                                   "q-value" = sum(object@qvalueObj$qvalues < x),
                                                   "local fdr" = sum(object@qvalueObj$lfdr < x)))
              colnames(counts) <- paste("<", cuts, sep = "")
              print(counts)
              cat("\n")
            }
          })
#' @rdname show
setMethod("show",
          signature = signature(object="deSet"),
          function(object) {
            cat('\n')
            cat('ExpressionSet Summary', '\n', '\n')
            tmp <- as(object, "ExpressionSet")
            print(tmp)
            cat('\n')
            cat('de Analysis Summary', '\n', '\n')
            cat('Total number of arrays:', ncol(exprs(object)), '\n')
            cat('Total number of probes:', nrow(exprs(object)), '\n', '\n')
            cat('Biological variables:', '\n')
            cat('\tNull Model: ')
            print(nullModel(object))
            cat('\tFull Model: ')
            print(fullModel(object))
            cat('\n')
            if (length(object@individual) != 0) {
              cat('Individuals:', '\n')
              ind <- as.numeric(object@individual)
              print(matrix(apply(((1:length(ind)) * t((ind))), 2, sum),
                           nrow=1))
              cat('\n')
            }
            cat('Expression data:', '\n')
            print(signif(exprs(object)[(1:min(5, nrow(exprs(object)))), ]),
                  digits = 3)
            cat('.......','\n','\n')
            if (!is.null(object@qvalueObj$pvalues)) {
              cuts <- c(0.0001, 0.001, 0.01, 0.025, 0.05, 0.10, 1)
              digits <- getOption("digits")
              cat("\nStatistical significance summary:\n")
              cat("pi0:", format(object@qvalueObj$pi0, digits = digits), "\n",
                  sep = "\t")
              cat("\n")
              cat("Cumulative number of significant calls:\n")
              cat("\n")
              counts <- sapply(cuts, function(x) c("p-value" = sum(object@qvalueObj$pvalues < x),
                                                   "q-value" = sum(object@qvalueObj$qvalues < x),
                                                   "local fdr" = sum(object@qvalueObj$lfdr < x)))
              colnames(counts) <- paste("<", cuts, sep="")
              print(counts)
              cat("\n")
            }
          })

#' @rdname apply_qvalue
setMethod("apply_qvalue",
          signature = signature(object="deSet"),
          function(object, ...) {
            if (length(object@qvalueObj) == 0) {
              stop("qvalueObj is empty- need to run either odp or lrt")
            }
            qvalueObj(object) <- qvalue(object@qvalueObj$pvalues, ...)
            validObject(object)
            object
          })

#' @rdname apply_sva
setMethod("apply_sva",
          signature = signature(object="deSet"),
          function(object, ...) {
            full.matrix <- object@full.matrix
            null.matrix <- object@null.matrix
            sv.sva <- sva(exprs(object),
                          mod0 = null.matrix,
                          mod = full.matrix, ...)$sv
            colnames(sv.sva) <- paste("SV", 1:ncol(sv.sva), sep="")
            pData(object) <- cbind(pData(object), sv.sva)
            fullModel(object) <- as.formula(paste("~",
                                                  paste(c(colnames(sv.sva),
                                                           attr(terms(fullModel(object)),
                                                                "term.labels")),
                                                         collapse=" + "),
                                                  sep=""))
            nullModel(object) <-  as.formula(paste("~",paste(c(colnames(sv.sva),
                                                                attr(terms(nullModel(object)),
                                                                     "term.labels")),
                                                             collapse=" + "),
                                                   sep=""))
            validObject(object)
            object
          })

#' @rdname apply_snm
setMethod("apply_snm",
          signature = signature(object="deSet"),
          function(object, int.var, ...) {
            full.matrix <- object@full.matrix
            null.matrix <- object@null.matrix
            full.matrix <- full.matrix - projMatrix(null.matrix) %*% full.matrix
            full.matrix <- as.matrix(rm.zero.cols(full.matrix))
            exprs(object) <- snm(exprs(object),
                                       bio.var = full.matrix,
                                       adj.var = null.matrix,
                                       int.var = int.var, ...)$norm.dat
            validObject(object)
            object
          })
