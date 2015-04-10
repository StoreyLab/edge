#' @rdname sType
setMethod("sType",
          signature = signature(object = "deFit"),
          function(object) {
            slot(object, "stat.type")
          })

#' @rdname betaCoef
setMethod("betaCoef",
          signature = signature(object = "deFit"),
          function(object) {
            slot(object, "beta.coef")
          })
#' @rdname resFull
setMethod("resFull",
          signature = signature(object = "deFit"),
          function(object) {
            slot(object, "res.full")
          })
#' @rdname resNull
setMethod("resNull",
          signature = signature(object = "deFit"),
          function(object) {
            slot(object, "res.null")
          })
#' @rdname fitFull
setMethod("fitFull",
          signature = signature(object = "deFit"),
          function(object) {
            slot(object, "fit.full")
          })
#' @rdname fitNull
setMethod("fitNull",
          signature = signature(object = "deFit"),
          function(object) {
            slot(object, "fit.null")
          })
#' @rdname fullModel
setMethod("fullModel",
          signature = signature(object = "deSet"),
          function(object) {
            slot(object, "full.model")
          })
#' @rdname nullModel
setMethod("nullModel",
          signature = signature(object = "deSet"),
          function(object) {
            slot(object, "null.model")
          })
#' @rdname fullMatrix
setMethod("fullMatrix",
          signature = signature(object = "deSet"),
          function(object) {
            slot(object, "full.matrix")
          })
#' @rdname nullMatrix
setMethod("nullMatrix",
          signature = signature(object = "deSet"),
          function(object) {
            slot(object, "null.matrix")
          })
#' @rdname individual
setMethod("individual",
          signature = signature(object = "deSet"),
          function(object) {
            slot(object, "individual")
          })
#' @rdname qvalueObj
setMethod("qvalueObj",
          signature = signature(object = "deSet"),
          function(object) {
            slot(object, "qvalueObj")
          })
