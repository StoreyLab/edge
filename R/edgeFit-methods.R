setMethod("summary",
          signature=signature(object="edgeFit"),
          function(object) {
            #SSE <- apply(object@res.full, 1, FUN=function(x) sum(x.^2))
            #SST <- object
            cat('\n'); cat('edgeFit Summary', '\n', '\n')
            cat('Models:', '\n')
            print(object@fitted.models)     
            cat('\nbeta.coef:', '\n')
            print(signif(object@beta.coef[(1:min(5, nrow(object@beta.coef))), ]), digits=3)
            cat('\nstat.type:', '\n')
            print(object@stat.type)
          })

setMethod("show",
          signature=signature(object="edgeFit"),
          function(object) {
            cat('\n'); cat('edgeFit Summary', '\n', '\n')
            cat('Models:', '\n')
            print(object@fitted.models)     
            cat('fit.full:', '\n')
            print(signif(object@fit.full[(1:min(2, nrow(object@fit.full))), ]), digits=3) 
            cat('\nfit.null:', '\n')
            print(signif(object@fit.null[(1:min(2, nrow(object@fit.null))), ]), digits=3) 
            cat('\nres.full:', '\n')
            print(signif(object@res.full[(1:min(2, nrow(object@res.full))), ]), digits=3) 
            cat('\nres.null:', '\n')
            print(signif(object@res.null[(1:min(2, nrow(object@res.null))), ]), digits=3)
            cat('\nbeta.coef:', '\n')
            print(signif(object@beta.coef[(1:min(5, nrow(object@beta.coef))), ]), digits=3)
            cat('\nstat.type:', '\n')
            print(object@stat.type)
          })