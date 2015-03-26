#' @rdname fullModel
setMethod("fullModel", 
          signature = signature(object = "edgeSet"), 
          function(object) {
            slot(object, "full.model")
          })
#' @rdname nullModel
setMethod("nullModel", 
          signature = signature(object = "edgeSet"), 
          function(object) {  
            slot(object, "null.model")
          })   
#' @rdname fullMatrix
setMethod("fullMatrix", 
          signature = signature(object = "edgeSet"), 
          function(object) {
            slot(object, "full.matrix")
          })
#' @rdname nullMatrix
setMethod("nullMatrix", 
          signature = signature(object = "edgeSet"), 
          function(object) {
            slot(object, "null.matrix")
          })   
#' @rdname individual
setMethod("individual", 
          signature = signature(object = "edgeSet"), 
          function(object) {
            slot(object, "individual")
          })     
#' @rdname qvalueObj
setMethod("qvalueObj", 
          signature = signature(object = "edgeSet"), 
          function(object) { 
            slot(object, "qvalueObj")
          })   
#' @rdname individual
setReplaceMethod("individual", 
                 signature = signature(object = "edgeSet"), 
                 function(object, value) {
                   object@individual <- value
                   validObject(object)
                   object            
                 })
#' @rdname qvalueObj
setReplaceMethod("qvalueObj",
                 signature = signature(object = "edgeSet"), 
                 function(object, value) {
                   object@qvalueObj <- value
                   validObject(object)
                   object            
                 })
#' @rdname fullModel
setReplaceMethod("fullModel", 
                 signature = signature(object = "edgeSet"), 
                 function(object, value) {
                   object@full.model <- value
                   fullMatrix(object) <- model.matrix(object = value, data = object)
                   validObject(object)
                   object
                 })
#' @rdname nullModel
setReplaceMethod("nullModel", 
                 signature = signature(object = "edgeSet"), 
                 function(object, value) {
                   object@null.model <- value
                   nullMatrix(object) <- model.matrix(object = value, data = object)
                   validObject(object)
                   object            
                 })
#' @rdname fullMatrix
setReplaceMethod("fullMatrix",
                 signature = signature(object = "edgeSet"), 
                 function(object, value) {
                   object@full.matrix <- value
                   validObject(object)
                   object            
                 })
#' @rdname nullMatrix
setReplaceMethod("nullMatrix", 
                 signature = signature(object = "edgeSet"), 
                 function(object, value) {
                   object@null.matrix <- value
                   validObject(object)
                   object            
                 })
#' @rdname modelFits
setMethod("modelFits", 
          signature = signature(object = "edgeFit"), 
          function(object) {  
            slot(object, "fitted.models")
          }) 
#' @rdname sType
setMethod("sType", 
          signature = signature(object = "edgeFit"), 
          function(object) {  
            slot(object, "stat.type")
          })   

#' @rdname betaCoef
setMethod("betaCoef", 
          signature = signature(object = "edgeFit"), 
          function(object) {  
            slot(object, "beta.coef")
          })
#' @rdname resFull
setMethod("resFull", 
          signature = signature(object = "edgeFit"), 
          function(object) {  
            slot(object, "res.full")
          }) 
#' @rdname resNull
setMethod("resNull", 
          signature = signature(object = "edgeFit"), 
          function(object) { 
            slot(object, "res.null")
          })
#' @rdname fitFull
setMethod("fitFull", 
          signature = signature(object = "edgeFit"), 
          function(object) {  
            slot(object, "fit.full")
          })  
#' @rdname fitNull
setMethod("fitNull", 
          signature = signature(object = "edgeFit"), 
          function(object) {  
            slot(object, "fit.null")
          })