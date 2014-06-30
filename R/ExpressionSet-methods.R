setAs("ExpressionSet", "edgeSet", function(from, to) updateOldExpSet(from, "edgeSet"))

updateOldExpSet <- function(from, toClass, ...) {  # to edgeSet
  # new object
  object <- new(toClass,
                assayData = from@assayData,
                phenoData = from@phenoData,
                featureData = annotatedDataFrameFrom(from@assayData,
                                                     byrow  = TRUE),
                experimentData = from@experimentData,
                annotation = from@annotation)
  validObject(object)
  object
}

setMethod("edgeSet", 
          signature = signature(object = "ExpressionSet"),
          function(object,
                   full.model,
                   null.model,
                   individual = NULL) {
            edgeObj <- as(object, "edgeSet")
            # Input checks
            if (!is.null(individual)) {
              if (length(individual) != ncol(exprs(object))) {
                stop("ind must be the same length as the number of arrays")
              }
            } 
            if (missing(full.model) || missing(null.model)) {
              stop("provide both full and null models")
            } 
            createSet(edgeObj,
                      nMod = null.model, 
                      fMod= full.model, 
                      ind = individual)
          })
