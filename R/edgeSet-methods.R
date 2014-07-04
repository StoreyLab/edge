setMethod("edgeFit",        
          "edgeSet",
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
              full.matrix <- full.matrix - H.null %*% full.matrix # subtract null model from full
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
            
            efObj <- new("edgeFit", fit.full = fitFull, fit.null = fitNull, 
                         dH.full = dHFull, res.full = resFull, res.null = resNull, 
                         beta.coef = B.coef, 
                         #fitted.models = list(full.model = object@full.model,
                        #                      null.model = object@null.model),
                         stat.type = stat.var)
            return(efObj) 
          })
setMethod("odp", 
          signature = signature(object = "edgeSet", obj.edgeFit = "missing"),
          function(object, obj.edgeFit, odp.parms = NULL, bs.its = 100, 
                   n.mods = 50, seed = NULL, verbose = TRUE, ...)  {
            obj.edgeFit <- edgeFit(object, 
                                   stat.type = "odp", ...)  
            results <- odp(object, obj.edgeFit, 
                           odp.parms = odp.parms,
                           n.mods = n.mods, 
                           bs.its = bs.its, 
                           seed = seed,
                           verbose = verbose, ...)
            return(results)
          })
setMethod("odp",
          signature = signature(object = "edgeSet", obj.edgeFit = "edgeFit"),
          function(object, obj.edgeFit, odp.parms = NULL, bs.its = 100, 
                   n.mods = 50, seed = NULL, verbose = TRUE, ...) {
            if (!is.null(seed)) {
              set.seed(seed)
            }     
            if (is.null(odp.parms)) {
              odp.parms <- klClust(object, 
                                   n.mods = n.mods, ...)
            } else if (sum(!(names(odp.parms) %in% c("mu.full", "sig.full", 
                                                     "mu.null", "sig.null", 
                                                     "n.per.mod", 
                                                     "clustMembers"))) != 0) {
              stop("Not a correct ODP parameter list. See klClust documentation")
            }
            odp.stat <- odpStat(n.res = obj.edgeFit@res.null, 
                                clustParms = odp.parms)
            null.stat <- bootstrap(object = object,
                                   obs.fit = obj.edgeFit,
                                   clustParms = odp.parms, 
                                   bs.its = bs.its,
                                   verbose = verbose)
            pval <- empPvals(stat = odp.stat,
                             stat0 = null.stat, ...)
            qvalue.obj(object) <- qvalue(p = pval, ...)
            return(object)
          })
setMethod("lrt", 
          signature = signature(object = "edgeSet", obj.edgeFit = "missing"),
          function(object, obj.edgeFit, nullDistn = c("normal", "bootstrap"), 
                   bs.its = 100, seed = NULL, ...) {
            edge.fit <- edgeFit(object,
                                stat.type = "lrt")  
            results <- lrt(object, 
                           obj.edgeFit = edge.fit, 
                           nullDistn = nullDistn, 
                           bs.its = bs.its,
                           seed = seed, 
                           verbose = verbose, ...)
            return(results)
          })
setMethod("lrt",
          signature = signature(object = "edgeSet", obj.edgeFit = "edgeFit"),
          function(object, obj.edgeFit, nullDistn = c("normal", "bootstrap"), 
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
            stat <- lrtStat(resNull = obj.edgeFit@res.null, 
                            resFull = obj.edgeFit@res.full)
            # If nullDistn is normal then return p-values from F-test else
            # return empirical p-values from qvalue package
            if (nullDistn == "normal") {
              df1 <- nFull - nNull
              df2 <- n - nFull
              stat <- stat * df2 / df1
              pval <- 1 - pf(stat,
                             df1 = df1, 
                             df2 = df2) 
              qvalue.obj(object) <- qvalue(p = pval, ...)
              return(object)
            } else {
              null.stat <- bootstrap(object = object, 
                                     obs.fit = obj.edgeFit, 
                                     bs.its = bs.its, 
                                     verbose = verbose)
              pval <- empPvals(stat = stat, 
                               stat0 = null.stat, ...)
              qvalue.obj(object) <- qvalue(pval, ...)
              return(object)
            }    
          })
setMethod("klClust",
          signature = signature(object = "edgeSet", obj.edgeFit = "missing"),
          function(object, obj.edgeFit,  n.mods = 50, ...)  {
            obj.edgeFit <- edgeFit(object, stat.type = "odp")  
            results <- klClust(object, obj.edgeFit, 
                               n.mods = n.mods, ...)
            return(results)
          })
setMethod("klClust", 
          signature = signature(object = "edgeSet", obj.edgeFit = "edgeFit"),
          function(object, obj.edgeFit, n.mods = 50, ...) {
            mod.member <- klmod(obj.edgeFit,
                                n.mods = n.mods)
            return(mod.parms(obj.edgeFit, 
                             clMembers = mod.member))
          })
setMethod("summary",
          signature = signature(object="edgeSet"),
          function(object, ...) {
            cat('\n')
            cat('ExpressionSet Summary', '\n', '\n')
            tmp <- as(object, "ExpressionSet")
            print(tmp)
            cat('\n')
            cat('edge Analysis Summary', '\n', '\n')
            cat('Total number of arrays:', ncol(exprs(object)), '\n')
            cat('Total number of probes:', nrow(exprs(object)), '\n', '\n')
            #  cat('Full model degrees of freedom:', df.model(object)$df.full, '\n')
            #  cat('Null model degrees of freedom:',df.model(object)$df.null, '\n', '\n')
            cat('Biological variables:', '\n')
            cat('\tNull Model:')
            print(nullModel(object))
            cat('\n\tFull Model:')
            print(fullModel(object))
            cat('\n') 
            if (length(object@individual) != 0) {
              cat('Individuals:', '\n')
              ind <- as.numeric(object@individual)
              print(matrix(apply(((1:length(ind)) * t((ind))), 2, sum), 
                           nrow = 1))
              cat('\n')
            }
            cat('.......', '\n', '\n')
            if (!is.null(object@qvalue.obj$pvalues)) {
              cuts <- c(0.0001, 0.001, 0.01, 0.025, 0.05, 0.10, 1)
              digits <- getOption("digits")
              cat("\nStatistical significance summary:\n")
              cat("pi0:", format(object@qvalue.obj$pi0, digits = digits), "\n", sep = "\t")
              cat("\n")
              cat("Cumulative number of significant calls:\n")
              cat("\n")
              counts <- sapply(cuts, function(x) c("p-value" = sum(object@qvalue.obj$pvalues < x), 
                                                   "q-value" = sum(object@qvalue.obj$qvalues < x), 
                                                   "local fdr" = sum(object@qvalue.obj$lfdr < x)))
              colnames(counts) <- paste("<", cuts, sep = "")
              print(counts)
              cat("\n")
            }
          })
setMethod("show",
          signature = signature(object="edgeSet"),
          function(object) {
            cat('\n')
            cat('ExpressionSet Summary', '\n', '\n')
            tmp <- as(object, "ExpressionSet")
            print(tmp)
            cat('\n')
            cat('edge Analysis Summary', '\n', '\n')
            cat('Total number of arrays:', ncol(exprs(object)), '\n')
            cat('Total number of probes:', nrow(exprs(object)), '\n', '\n')
           # cat('Full model degrees of freedom:', dfModel(object)$df.full, '\n')
           # cat('Null model degrees of freedom:',dfModel(object)$df.null, '\n', '\n')
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
            print(signif(exprs(object)[(1:min(5, nrow(exprs(object)))), ]), digits = 3) 
            cat('.......','\n','\n')
            if (!is.null(object@qvalue.obj$pvalues)) {
              cuts <- c(0.0001, 0.001, 0.01, 0.025, 0.05, 0.10, 1)
              digits <- getOption("digits")
              cat("\nStatistical significance summary:\n")
              cat("pi0:", format(object@qvalue.obj$pi0, digits = digits), "\n", sep = "\t")
              cat("\n")
              cat("Cumulative number of significant calls:\n")
              cat("\n")
              counts <- sapply(cuts, function(x) c("p-value" = sum(object@qvalue.obj$pvalues < x), 
                                                   "q-value" = sum(object@qvalue.obj$qvalues < x), 
                                                   "local fdr" = sum(object@qvalue.obj$lfdr < x)))
              colnames(counts) <- paste("<", cuts, sep="")
              print(counts)
              cat("\n")
            }
          })
setMethod("edgeQvalue",
          signature = signature(object="edgeSet"),
          function(object, ...) {
            if (length(object@qvalue.obj) == 0) {
              stop("qvalue.obj is empty- need to run either odp or lrt")
            }
            qvalue.obj(object) <- qvalue(object@qvalue.obj$pvalues, ...)
            validObject(object)
            object
          })

setMethod("edgeSVA",
          signature = signature(object="edgeSet"),
          function(object, ...) {
            full.matrix <- object@full.matrix
            null.matrix <- object@null.matrix
          #  full.matrix <- full.matrix - projMatrix(null.matrix) %*% full.matrix 
          #  full.matrix <- as.matrix(rm.zero.cols(full.matrix))
            sv.sva <- sva(exprs(object),
                          mod0 = null.matrix,
                          mod = full.matrix, ...)$sv
            nullMatrix(object) <- cbind(sv.sva, object@null.matrix)
            fullMatrix(object) <- cbind(sv.sva, object@full.matrix)
            validObject(object)
            object
          })

setMethod("edgeSNM",
          signature=signature(object="edgeSet"),
          function(object, int.var, ...) {
            full.matrix <- object@full.matrix
            null.matrix <- object@null.matrix
            full.matrix <- full.matrix - projMatrix(null.matrix) %*% full.matrix
            full.matrix <- as.matrix(rm.zero.cols(full.matrix))
            nullMatrix(object) <- snm(exprs(object),
                                       bio.var = full.matrix,
                                       adj.var = null.matrix, 
                                       int.var = int.var, ...)$norm.dat
            validObject(object)
            object
          })