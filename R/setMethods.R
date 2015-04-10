#' @rdname individual
setReplaceMethod("individual",
                 signature = signature(object = "deSet"),
                 function(object, value) {
                   object@individual <- value
                   validObject(object)
                   object
                 })
#' @rdname qvalueObj
setReplaceMethod("qvalueObj",
                 signature = signature(object = "deSet"),
                 function(object, value) {
                   object@qvalueObj <- value
                   validObject(object)
                   object
                 })
#' @rdname fullModel
setReplaceMethod("fullModel",
                 signature = signature(object = "deSet"),
                 function(object, value) {
                   object@full.model <- value
                   fullMatrix(object) <- model.matrix(object = value, data = object)
                   validObject(object)
                   object
                 })
#' @rdname nullModel
setReplaceMethod("nullModel",
                 signature = signature(object = "deSet"),
                 function(object, value) {
                   object@null.model <- value
                   nullMatrix(object) <- model.matrix(object = value, data = object)
                   validObject(object)
                   object
                 })
#' @rdname fullMatrix
setReplaceMethod("fullMatrix",
                 signature = signature(object = "deSet"),
                 function(object, value) {
                   object@full.matrix <- value
                   validObject(object)
                   object
                 })
#' @rdname nullMatrix
setReplaceMethod("nullMatrix",
                 signature = signature(object = "deSet"),
                 function(object, value) {
                   object@null.matrix <- value
                   validObject(object)
                   object
                 })
