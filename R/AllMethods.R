setMethod("fullModel", 
          signature = signature(object = "edgeSet"), 
          function(object) {
            slot(object, "full.model")
          })
setMethod("nullModel", 
          signature = signature(object = "edgeSet"), 
          function(object) {  
            slot(object, "null.model")
          })   
setMethod("fullMatrix", 
          signature = signature(object = "edgeSet"), 
          function(object) {
            slot(object, "full.matrix")
          })
setMethod("nullMatrix", 
          signature = signature(object = "edgeSet"), 
          function(object) {
            slot(object, "null.matrix")
          })   

setMethod("individual", 
          signature = signature(object = "edgeSet"), 
          function(object) {
            slot(object, "individual")
          })     
setMethod("qvalueObj", 
          signature = signature(object = "edgeSet"), 
          function(object) { 
            slot(object, "qvalueObj")
          })   
setReplaceMethod("individual", 
                 signature = signature(object = "edgeSet"), 
                 function(object, value) {
                   object@individual <- value
                   validObject(object)
                   object            
                 })
setReplaceMethod("qvalueObj",
                 signature = signature(object = "edgeSet"), 
                 function(object, value) {
                   object@qvalueObj <- value
                   validObject(object)
                   object            
                 })

setReplaceMethod("fullModel", 
                 signature = signature(object = "edgeSet"), 
                 function(object, value) {
                   object@full.model <- value
                   validObject(object)
                   object
                 })

setReplaceMethod("nullModel", 
                 signature = signature(object = "edgeSet"), 
                 function(object, value) {
                   object@null.model <- value
                   validObject(object)
                   object            
                 })
setReplaceMethod("fullMatrix",
                 signature = signature(object = "edgeSet"), 
                 function(object, value) {
                   object@full.matrix <- value
                   validObject(object)
                   object            
                 })
setReplaceMethod("nullMatrix", 
                 signature = signature(object = "edgeSet"), 
                 function(object, value) {
                   object@null.matrix <- value
                   validObject(object)
                   object            
                 })

# edgeFit get-methods
setMethod("modelFits", 
          signature = signature(object = "edgeFit"), 
          function(object) {  
            slot(object, "fitted.models")
          }) 
setMethod("sType", 
          signature = signature(object = "edgeFit"), 
          function(object) {  
            slot(object, "stat.type")
          })   
setMethod("betaCoef", 
          signature = signature(object = "edgeFit"), 
          function(object) {  
            slot(object, "beta.coef")
          })
setMethod("resFull", 
          signature = signature(object = "edgeFit"), 
          function(object) {  
            slot(object, "res.full")
          })  
setMethod("resNull", 
          signature = signature(object = "edgeFit"), 
          function(object) { 
            slot(object, "res.null")
          })
setMethod("fitFull", 
          signature = signature(object = "edgeFit"), 
          function(object) {  
            slot(object, "fit.full")
          })  
setMethod("fitNull", 
          signature = signature(object = "edgeFit"), 
          function(object) {  
            slot(object, "fit.null")
          })