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