#' Get model equations from edgeFit object
#'
#' Performs F-test (Generalized Likelihood Ratio)
#'
#' \code{lrt} performs Generalized Likelihood Ratio on the full and null models to determine p-values (either from a F-distribution or
#' from simulated null test statistics). If the null distribution is calculated using "bootstrap" then
#' residuals from the alternative model are resampled and added to the null model to simulate
#' a distribution where there is no differential expression. Otherwise, the default
#' input is "normal" and the assumption is that the dataset follows an F-distribution. Note that
#' if "normal" is chosen then the sum of squares should be statistically independent
#' and genes should be normally distributed with the same variance. The \code{qvalue} function
#' is used to retrieve the \code{qvalue} object. See \code{\link{qvalue}} for more information.
#'
#' @param object \code{\linkS4class{edgeSet}}
#' @param obj.edgeFit edgeFit object- S4 class
#' @param nullDistn character- either "normal" or "bootstrap", If "normal",
#' use p-values from ANOVA else if "bootstrap" p-values will
#' be determined from function \code{\link{empPvals}}.
#' Default is "normal".
#' @param bs.its numeric- number of null statistics generated (only applicable
#' for "bootstrap" method. Default is 100.
#' @param seed numeric- set the seed value.
#' @param verbose boolean- print iterations for bootstrap method. Default is TRUE.
#' @param ... Additional arguments for \code{qvalue} and \code{empPvals} function.
#'
#' @note Fits a linear regression to each gene from function
#' \code{\link{edgeFit}} and then performs an ANOVA F-test on the null and full
#' model fits. If nullDistn "normal" is chosen then the F-distribution is used to determine
#' p-values else if "bootstrap" is chosen then empirical p-values will be
#' determined from \code{\link{qvalue}} package (see \code{\link{empPvals}}).
#' To retrieve linear regression results, see \code{\link{edgeFit}}.
#'
#' @author John Storey, Andrew Bass
#'
#' @return \code{lrt} returns a \code{\linkS4class{edgeSet}} object with slot qvalue.obj as \code{\link{qvalue}} object.
#'
#' @examples
#' # Create ExpressionSet object from kidney dataset
#' library(splines)
#' data(kidney)
#' sex <- kidney$sex[kidney$tissue == "c"]
#' age <- kidney$age[kidney$tissue == "c"]
#' kidexpr <- log(kidney$kidexpr[, kidney$tissue == "c"] + 10)
#' expSet <- ExpressionSet(assayData = kidexpr, phenoData = as(data.frame(sex = sex, age = age), "AnnotatedDataFrame"))
#'
#' # Create Models
#' nModel <- ~sex
#' fModel <- ~sex+ ns(age, df=3, intercept=FALSE)
#'
#' # Create edgeSet
#' edgeObj <- edgeSet(expSet, full.model=fModel, null.model=nModel)
#'
#' # lrt method
#' edge.lrt <- lrt(edgeObj, nullDistn="normal")
#'
#' # To generate pvalues from bootstrap
#' edge.lrt.bs <- lrt(edgeObj, nullDistn="bootstrap")
#'
#' # Input an edgeFit object but not necessary
#' edge.fit <- edgeFit(edgeObj, stat.type="lrt")
#' edge.lrt.ef <- lrt(edgeObj, obj.edgeFit=edge.fit)
#'
#' @references
#' Storey JD, Xiao W, Leek JT, Tompkins RG, and Davis RW. (2005) Significance analysis of time course microarray experiments. Proceedings of the National Academy of Sciences, 102: 12837-12842.
#'
#' @seealso \code{\link{edgeSet}}, \code{\link{odp}}, \code{\link{edgeFit}}
#'
#' @keywords lrt
#'
#' @rdname lrt-method
#' @import Biobase
#' @import methods
#' @exportMethod lrt
setGeneric("lrt", function(object, obj.edgeFit,
                           nullDistn=c("normal","bootstrap"), bs.its=100,
                           seed=NULL, verbose=TRUE, ...)
  standardGeneric("lrt"))


#' Performs Optimal Discovery Procedure (ODP) on edgeSet object
#'
#' \code{odp} performs the Optimal Discovery Procedure, which is a new approach for optimally performing
#' many hypothesis tests in a high-dimensional study. The method has been introduced recently (Storey 2007)
#' and has been theoretically shown to optimally perform multiple significance tests. When testing a feature,
#' information from all the features is utilized when testing for significance of a feature. It guarentees to maximize the
#' number of expected true positive results for each fixed number of expected false positive results which is related to FDR. An
#' ODP score/statistic is determined for each gene (refer to Storey 2007 for more details) and the null
#' distribution is generated through a resampling method called bootstrap where residuals from the alternative
#' model are resampled and added back to the null model to simulate the case where there is no differential
#' expression.
#'
#' @param object \code{\linkS4class{edgeSet}}
#' @param obj.edgeFit edgeFit object.
#' @param odp.parms list- clustering parameters.
#' @param bs.its numeric- number of null statistics generated. Default
#' is 100.
#' @param n.mods numeric- number of clusters.
#' @param seed numeric- set the seed value.
#' @param verbose boolean- print iterations for bootstrap method. Default is TRUE.
#' @param ... Additional arguments for \code{klClust}, \code{qvalue} and \code{empPvals}.
#'
#'
#' @details The full ODP estimator computationally grows quadratically with respect to
#' the number of genes. This is because the likelihood function of each gene be
#' calculated across all parameter estimates of the rest of the genes. This becomes
#' computationally infeasible at a certain point.
#' Therefore, an alternative method called mODP is used which
#' has been shown to provide results that are very similar. mODP utilizes a k-means
#' clustering algorithm where genes are assigned to a cluster based on the Kullback-Leiber distance.
#' Each gene is assigned an module-average parameter to calculate the ODP score
#' and it reduces the quadratic time to linear (See Woo, Leek and Storey 2010).
#'
#' @note The null and full models are fitted using a linear regression.
#' To retrieve linear regression results, see \code{\link{edgeFit}}. The p-values
#' are determined by function \code{\link{empPvals}}.
#'
#' @return \code{odp} returns an \code{\linkS4class{edgeSet}} object with slot qvalue.obj as
#' \code{\link{qvalue}} object.
#'
#' @examples
#' # Create ExpressionSet object from kidney dataset
#' library(splines)
#' data(kidney)
#' sex <- kidney$sex[kidney$tissue == "c"]
#' age <- kidney$age[kidney$tissue == "c"]
#' kidexpr <- log(kidney$kidexpr[, kidney$tissue == "c"] + 10)
#' expSet <- ExpressionSet(assayData = kidexpr, phenoData = as(data.frame(sex = sex, age = age), "AnnotatedDataFrame"))
#'
#' # Create Models
#' nModel <- ~sex
#' fModel <- ~sex+ ns(age, df=3, intercept=FALSE)
#'
#' # Create edgeSet
#' edgeObj <- edgeSet(expSet, full.model=fModel, null.model=nModel)
#'
#' # ODP method
#' edge.odp <- odp(edgeObj)
#'
#' # Change the number of clusters
#' edge.lrt.bs <- lrt(edgeObj, bs.its=10)
#'
#' # Input an edgeFit object or ODP parameters... not necessary
#' edge.fit <- edgeFit(edgeObj, stat.type="odp")
#' edge.parms <- klClust(edgeObj)
#' edge.odp.ef <- odp(edgeObj, obj.edgeFit=edge.fit, odp.parms=edge.parms)
#'
#' @references
#' Storey JD. (2007) The optimal discovery procedure: A new approach to simultaneous significance testing. Journal of the Royal Statistical Society, Series B, 69: 347-368.
#'
#' Storey JD, Dai JY, and Leek JT. (2007) The optimal discovery procedure for large-scale significance testing, with applications to comparative microarray experiments. Biostatistics, 8: 414-432.
#'
#' Woo S, Leek JT, Storey JD (2010) A computationally efficient modular optimal discovery procedure. Bioinformatics, 27(4): 509-515.
#'
#' @author John Storey, Andrew Bass
#'
#' @seealso \code{\link{klClust}}, \code{\link{edgeSet}} and \code{\link{edgeFit}}
#'
#' @keywords odp
#'
#' @rdname odp-method
#' @useDynLib edge odpScoreCluster kldistance
#' @import qvalue MASS splines
#' @exportMethod odp
setGeneric("odp", function(object, obj.edgeFit, odp.parms=NULL, bs.its=100,
                           n.mods=50, seed=NULL, verbose=TRUE, ...)
  standardGeneric("odp"))


#' Clustering (mODP) method on edgeSet object
#'
#' \code{klClust} is an implementation of mODP that assigns genes to modules based
#' off of the Kullback-Leibler distance.
#'
#' @param object \code{\linkS4class{edgeSet}}
#' @param obj.edgeFit edgeFit object.
#' @param n.mods numeric- number of clusters.
#' @param \dots additional parameters
#'
#'
#' @details mODP utilizes a k-means clustering algorithm where genes are assigned
#' to a cluster based on the Kullback-Leiber distance. Each gene is assigned an
#' module-average parameter to calculate the ODP score (See Woo, Leek and Storey 2010 for more details).
#' The mODP and full ODP produce near exact results but mODP has the advantage of being computationally
#' feasible.
#'
#' @note The results are generally insensitive to the number of modules after a
#' certain threshold of about K>=50.
#'
#' @return
#' \code{klClust} returns a list:
#' \itemize{
#' \item {mu.full: mean of clusters from full model}
#' \item {mu.null: mean of clusters from null model}
#' \item {sig.full: sd of clusters from full model}
#' \item {sig.null: sd of clusters from null model}
#' \item {n.per.mod: total members in each cluster}
#' \item {clustMembers: members for each gene}
#'
#' }
#'
#' @examples
#' # Create ExpressionSet object from kidney dataset
#' library(splines)
#' data(kidney)
#' sex <- kidney$sex[kidney$tissue == "c"]
#' age <- kidney$age[kidney$tissue == "c"]
#' kidexpr <- log(kidney$kidexpr[, kidney$tissue == "c"] + 10)
#' expSet <- ExpressionSet(assayData = kidexpr, phenoData = as(data.frame(sex = sex, age = age), "AnnotatedDataFrame"))
#'
#' # Create Models
#' nModel <- ~sex
#' fModel <- ~sex+ ns(age, df=3, intercept=FALSE)
#'
#' # Create edgeSet
#' edgeObj <- edgeSet(expSet, full.model=fModel, null.model=nModel)
#'
#' # ODP method
#' edge.clust <- klClust(edgeObj)
#'
#' # Change the number of clusters
#' edge.clust <- klClust(edgeObj, n.mods=10)
#'
#' # Input an edgeFit object or ODP parameters... not necessary
#' edge.fit <- edgeFit(edgeObj, stat.type="odp")
#' edge.clust.ef <- klClust(edgeObj, obj.edgeFit=edge.fit)
#'
#' @references
#' Storey JD. (2007) The optimal discovery procedure: A new approach to simultaneous significance testing. Journal of the Royal Statistical Society, Series B, 69: 347-368.
#'
#' Storey JD, Dai JY, and Leek JT. (2007) The optimal discovery procedure for large-scale significance testing, with applications to comparative microarray experiments. Biostatistics, 8: 414-432.
#'
#' Woo S, Leek JT, Storey JD (2010) A computationally efficient modular optimal discovery procedure. Bioinformatics, 27(4): 509-515.
#'
#' @author John Storey, Andrew Bass
#'
#' @seealso \code{\link{odp}}, \code{\link{edgeSet}} and \code{\link{edgeFit}}
#'
#' @rdname klClust-method
#' @keywords klClust
#' @exportMethod klClust
setGeneric("klClust", function(object, obj.edgeFit=NULL, n.mods=50, ...)
  standardGeneric("klClust"))

#' Linear Regression for full and null models
#'
#' \code{edgeFit} fits a linear model to each gene. Model fits can be either statistic type
#' odp (Optimal Discovery Procedure) or lrt (Likelihood Ratio Test).
#'
#' @param object \code{\linkS4class{edgeSet}} object.
#' @param stat.type character- type of statistic to be used. Either "lrt" or "odp".
#' Default is "lrt".
#'
#' @details If individual factors exists from \code{edgeSet} object, they are removed from the model to remove
#' the effects from individuals. If ODP method is implemented then the null model
#' is removed from the full model (see Storey 2007). Statistics of interest for functions \code{lrt}, \code{odp} and
#' \code{klClust} are the elements of the projection matrix, the fitted
#' expression data, the regression coefficients and the residuals of the expression data for both the null
#' and full models.
#'
#' @note \code{edgeFit} does not have to be called by the user to use
#' \code{odp}, \code{klClust} or \code{lrt} as it is an optional input and is
#' implemented in the methods if not specified. edgeFit object can be created
#' by the user if a different statistical implementation is required.
#'
#' @return \code{edgeFit} returns an \code{\linkS4class{edgeFit}} object.
#'
#' @examples
#' # Create ExpressionSet object from kidney dataset
#' library(splines)
#' data(kidney)
#' sex <- kidney$sex[kidney$tissue == "c"]
#' age <- kidney$age[kidney$tissue == "c"]
#' kidexpr <- log(kidney$kidexpr[, kidney$tissue == "c"] + 10)
#' expSet <- ExpressionSet(assayData = kidexpr, phenoData = as(data.frame(sex = sex, age = age), "AnnotatedDataFrame"))
#'
#' # Create Models
#' nModel <- ~sex
#' fModel <- ~sex + ns(age, df=3, intercept=FALSE)
#'
#' # edgeSet object
#' edgeObj <- edgeSet(expSet, full.model=fModel, null.model=nModel)
#'
#' # Retrieve statistics from linear regression for each gene
#' ef.lrt <- edgeFit(edgeObj, stat.type="lrt") # lrt method
#' ef.odp <- edgeFit(edgeObj, stat.type="odp") # odp method
#'
#' @references
#' Storey JD. (2007) The optimal discovery procedure: A new approach to simultaneous significance testing. Journal of the Royal Statistical Society, Series B, 69: 347-368.
#'
#' Storey JD, Dai JY, and Leek JT. (2007) The optimal discovery procedure for large-scale significance testing, with applications to comparative microarray experiments. Biostatistics, 8: 414-432.
#'
#' Storey JD, Xiao W, Leek JT, Tompkins RG, and Davis RW. (2005) Significance analysis of time course microarray experiments. Proceedings of the National Academy of Sciences, 102: 12837-12842.
#'
#' @seealso \code{\linkS4class{edgeFit}}, \code{\link{odp}} and
#' \code{\link{lrt}}
#'
#' @author John Storey, Andrew Bass
#' @rdname edgeFit-method
#' @keywords edgeFit
#' @exportMethod edgeFit
setGeneric("edgeFit",
           function(object, stat.type=c("lrt", "odp")) {
             standardGeneric("edgeFit")
           })

#' Creates object of class edgeSet
#'
#' Creates an edgeSet object that extends ExpressionSet.
#'
#' @param object \code{\link{ExpressionSet}} object
#' @param full.model formula- full model for experiment data.
#' @param null.model formula- null model for experiment data.
#' @param individual factor- information on individuals in experiment.
#'
#' @note It is essential that the null and full models have the same variables
#' as the ExpressionSet phenoType column names.
#'
#' @return \code{edgeSet} returns an \code{\linkS4class{edgeSet}} object.
#'
#' @examples
#' # Create ExpressionSet object from kidney dataset
#' library(splines)
#' data(kidney)
#' sex <- kidney$sex[kidney$tissue == "c"]
#' age <- kidney$age[kidney$tissue == "c"]
#' kidexpr <- log(kidney$kidexpr[, kidney$tissue == "c"] + 10)
#' expSet <- ExpressionSet(assayData = kidexpr, phenoData = as(data.frame(sex = sex, age = age), "AnnotatedDataFrame"))
#'
#' # Create Models
#' nModel <- ~sex
#' fModel <- ~sex+ ns(age, df=3, intercept=FALSE)
#'
#' # Create edgeSet
#' edgeObj <- edgeSet(expSet, full.model=fModel, null.model=nModel)
#'
#' # Optionally add individuals to experiment for kidney data not applicable
#' edgeObj <- edgeSet(expSet, full.model=fModel, null.model=nModel, ind=factor(1:72))
#'
#'
#' @seealso \code{\linkS4class{edgeSet}}, \code{\link{odp}} and \code{\link{lrt}}
#'
#' @author John Storey, Andrew Bass
#' @rdname edgeSet-method
#' @keywords edgeSet
#'
#' @exportMethod edgeSet
setGeneric("edgeSet", function(object, full.model, null.model, individual=NULL) standardGeneric("edgeSet"))

#' Estimate the q-values for a given set of p-values in edgeSet object
#'
#' Runs \code{qvalue} on edgeSet object based on p-values determined from functions
#' \code{lrt} or \code{odp}
#'
#' @param object \code{\linkS4class{edgeSet}} object
#' @param ... Additional arguments for \code{\link{qvalue}}
#'
#' @return \code{edgeQvalue} returns an \code{\linkS4class{edgeSet}} object.
#'
#' @examples
#' # Create ExpressionSet object from kidney dataset
#' library(splines)
#' data(kidney)
#' sex <- kidney$sex[kidney$tissue == "c"]
#' age <- kidney$age[kidney$tissue == "c"]
#' kidexpr <- log(kidney$kidexpr[, kidney$tissue == "c"] + 10)
#' expSet <- ExpressionSet(assayData = kidexpr, phenoData = as(data.frame(sex = sex, age = age), "AnnotatedDataFrame"))
#'
#' # Create Models
#' nModel <- ~sex
#' fModel <- ~sex+ ns(age, df=3, intercept=FALSE)
#'
#' # Create edgeSet
#' edgeObj <- edgeSet(expSet, full.model=fModel, null.model=nModel)
#'
#' # Run lrt (or odp) and edgeQvalue
#' edge.lrt <- lrt(edgeObj)
#' edge.lrt.new <- edgeQvalue(edge.lrt, fdr.level=0.05, pi0.method="bootstrap", adj=1.2)
#'
#' @rdname edgeQvalue-method
#'
#' @references
#' Storey JD and Tibshirani R. (2003) Statistical significance for genome-wide studies. Proceedings of the National Academy of Sciences, 100: 9440-9445
#'
#' @seealso \code{\linkS4class{edgeSet}}, \code{\link{odp}} and \code{\link{lrt}}
#'
#' @author John Storey, Andrew Bass
#'
#' @keywords edgeQvalue
#' @rdname edgeQvalue-method
#' @exportMethod edgeQvalue
setGeneric("edgeQvalue", function(object, ...)
  standardGeneric("edgeQvalue"))

#' Estimate surrogate variables 
#'
#' Runs \code{sva} on edgeSet object based on the null and full models in \code{\linkS4class{edgeSet}}
#'
#' @param object \code{\linkS4class{edgeSet}} object
#' @param ... Additional arguments for \code{\link{sva}}
#'
#' @return \code{edgeSVA} returns an \code{\linkS4class{edgeSet}} object.
#' 
#' @examples
#' # Create ExpressionSet object from kidney dataset
#' library(splines)
#' data(kidney)
#' sex <- kidney$sex[kidney$tissue == "c"]
#' age <- kidney$age[kidney$tissue == "c"]
#' kidexpr <- log(kidney$kidexpr[, kidney$tissue == "c"] + 10)
#' expSet <- ExpressionSet(assayData = kidexpr, phenoData = as(data.frame(sex = sex, age = age), "AnnotatedDataFrame"))
#'
#' # Create Models
#' nModel <- ~sex
#' fModel <- ~sex+ ns(age, df=3, intercept=FALSE)
#'
#' # Create edgeSet
#' edgeObj <- edgeSet(expSet, full.model=fModel, null.model=nModel)
#' edgeObj.sva <- edgeSVA(edgeObj)
#'
#' @seealso \code{\linkS4class{edgeSet}}, \code{\link{odp}} and \code{\link{lrt}}
#'
#' @references
#' Leek JT, Storey JD (2007) Capturing Heterogeneity in Gene Expression Studies by Surrogate Variable Analysis. PLoS Genet 3(9): e161. doi:10.1371/journal.pgen.0030161
#'
#' @author John Storey, Andrew Bass
#' @import sva
#' @keywords edgeSVA
#' @rdname edgeSVA-method
#' @exportMethod edgeSVA
setGeneric("edgeSVA", function(object, ...)
  standardGeneric("edgeSVA"))

#' Supervised normalization of data in edgeSet object
#'
#' Runs \code{snm} on edgeSet object based on the null and full models in \code{\linkS4class{edgeSet}}
#'
#' @param object \code{\linkS4class{edgeSet}} object
#' @param int.var data frame- intensity-dependent effects.
#' @param ... Additional arguments for \code{\link{snm}}
#'
#' @return \code{edgeSNM} returns an \code{\linkS4class{edgeSet}} object.
#'
#' @references
#' Mechan BH, Nelson PS, Storey JD. Supervised normalization of microarrays. Bioinformatics 2010;26:1308-15
#'
#' @examples
#' # Simulate data
#' require(snm)
#' singleChannel <- sim.singleChannel(12345)
#' # Create edgeSet object using edgeModel (can use ExpressionSet see manual)
#' edgeObj <- edgeModel(data=singleChannel$raw.data, sampling="static", bio.var=singleChannel$bio.var, adj.var=singleChannel$adj.var)
#'
#'
#' # Run SNM using intensity-dependent adjustment variable
#' nEdgeObj <- edgeSNM(edgeObj, int.var=singleChannel$int.var, verbose=FALSE, num.iter=1)
#'
#' @seealso \code{\linkS4class{edgeSet}}, \code{\link{odp}} and \code{\link{lrt}}
#'
#' @author John Storey, Andrew Bass
#' @rdname edgeSNM-method
#' @keywords edgeSNM
#' @import snm
#' @exportMethod edgeSNM
setGeneric("edgeSNM", function(object, int.var, ...) standardGeneric("edgeSNM"))
#' Full model equation from edgeSet object
#' 
#' These generic functions access and set the full model for 
#' \code{\linkS4class{edgeSet}} object.
#'
#' @param object \code{\linkS4class{edgeSet}}
#' @param value \code{formula}: The experiment design for the full model.
#' 
#' @usage fullModel(object)
#' 
#' @return \code{fullModel} returns the formula for the full model.
#'    
#' @author John Storey, Andrew Bass
#' 
#' @seealso \code{\linkS4class{edgeSet}}
#' 
#' @keywords fullModel, fullModel<-
#' 
#' @rdname fullModel-method
#' 
#' @exportMethod fullModel
setGeneric("fullModel", function(object) standardGeneric("fullModel"))

#' @rdname fullModel-method
#' 
#' @exportMethod fullModel<-
setGeneric("fullModel<-", function(object, value) {
  standardGeneric("fullModel<-") 
})

#' Null model equation from edgeSet object
#'
#' These generic functions access and set the null model for 
#' \code{\linkS4class{edgeSet}} object.
#'
#' @param object \code{\linkS4class{edgeSet}}
#' @param value \code{formula}: The experiment design for the null model.
#' 
#' @usage nullModel(object)
#' 
#' @return \code{nullModel} returns the formula for the null model.
#'
#' @author John Storey, Andrew Bass
#' 
#' @seealso \code{\linkS4class{edgeSet}}
#' 
#' @keywords nullModel, nullModel<-
#' 
#' @rdname nullModel-method
#' @exportMethod nullModel
setGeneric("nullModel", function(object) standardGeneric("nullModel"))

#' @rdname nullModel-method
#' 
#' @exportMethod nullModel<-
setGeneric("nullModel<-", function(object, value) {
  standardGeneric("nullModel<-") 
})

#' Matrix representation of null model
#'
#' These generic functions access and set the null matrix for 
#' \code{\linkS4class{edgeSet}} object.
#'
#' @param object \code{\linkS4class{edgeSet}}
#' @param value \code{matrix}: null model matrix where columns are covariates and rows are observations
#' 
#' @usage nullMatrix(object)
#' 
#' @return \code{nullMatrix} returns the value of the null model matrix.
#'
#' @author John Storey, Andrew Bass
#' 
#' @seealso \code{\linkS4class{edgeSet}}, \code{\link{fullModel}} and 
#' \code{\link{fullModel}}
#' 
#' @keywords nullMatrix, nullMatrix<-
#' 
#' @rdname nullMatrix-method   
#' 
#' @exportMethod nullMatrix
setGeneric("nullMatrix", function(object) standardGeneric("nullMatrix"))

#' @rdname nullMatrix-method
#'  
#' @exportMethod nullMatrix<-
setGeneric("nullMatrix<-", function(object, value) {
  standardGeneric("nullMatrix<-") 
})

#' Matrix representation of full model
#'
#' These generic functions access and set the full matrix for 
#' \code{\linkS4class{edgeSet}} object.
#'
#' @param object \code{\linkS4class{edgeSet}}
#' @param value \code{matrix}: full model matrix where the columns are the covariates amd rows are observations
#' 
#' @usage fullMatrix(object)
#' 
#' @return \code{fullMatrix} returns the value of the full model matrix.
#'
#' @author Andrew Bass, John Storey
#' 
#' @seealso \code{\linkS4class{edgeSet}}, \code{\link{fullModel}}
#' 
#' @keywords fullMatrix, fullMatrix<-
#' 
#' @rdname fullMatrix-method
#' 
#' @exportMethod fullMatrix
setGeneric("fullMatrix", function(object) standardGeneric("fullMatrix"))

#' @rdname fullMatrix-method
#' 
#' @exportMethod fullMatrix<-
setGeneric("fullMatrix<-", function(object, value) {
  standardGeneric("fullMatrix<-") 
})


#' qvalue object in experiment
#'
#' These generic functions access and set the \code{qvalue} object in the
#' \code{\linkS4class{edgeSet}} object.
#'
#' @usage qvalue.obj(object)
#' 
#' @param object \code{\linkS4class{edgeSet}}
#' @param value S3 \code{object}: \code{\link{qvalue}} object
#' 
#' @return  \code{qvalue.obj} returns a \code{\link{qvalue}} object.
#'     
#' @author John Storey, Andrew Bass
#' 
#' @seealso \code{\link{lrt}}, \code{\link{odp}} and 
#' \code{\linkS4class{edgeSet}}
#' 
#' @keywords qvalue.obj, qvalue.obj<-
#' 
#' @rdname qvalue.obj-method  
#' 
#' @exportMethod qvalue.obj                
setGeneric("qvalue.obj", function(object) standardGeneric("qvalue.obj"))

#' @rdname qvalue.obj-method
#' 
#' @exportMethod qvalue.obj<-
setGeneric("qvalue.obj<-", function(object, value) {
  standardGeneric("qvalue.obj<-") 
})                      

#' Individuals utilized in experiment
#'
#' These generic functions access and set the individual slot in 
#' \code{\linkS4class{edgeSet}}.
#'
#' @usage individual(object)
#' 
#' @param object \code{\linkS4class{edgeSet}}
#' @param value \code{factor}: individual identifier for each observation
#' 
#' @return \code{individual} returns information regarding individuals 
#' in the experiment.
#'    
#' @author John Storey, Andrew Bass
#' 
#' @seealso \code{\linkS4class{edgeSet}}
#' 
#' @keywords individual, individual<-
#' 
#' @rdname individual-method
#' 
#' @exportMethod individual
setGeneric("individual", function(object) standardGeneric("individual"))

#' @rdname individual-method
#' 
#' @exportMethod individual<-
setGeneric("individual<-", function(object, value) {
  standardGeneric("individual<-") 
})  
#' Fitted models for an \code{\linkS4class{edgeFit}} object.
#'
#' @param object \code{\linkS4class{edgeFit}} object
#' 
#' @usage modelFits(object)
#' 
#' @return \code{modelFits} returns the formula for the fitted models.
#'
#' @author John Storey, Andrew Bass
#' 
#' @seealso \code{\link{edgeFit}}, \code{\linkS4class{edgeSet}}
#' 
#' @keywords modelFits
#' 
#' @rdname modelFits-method
#' 
#' @exportMethod modelFits
setGeneric("modelFits", function(object) standardGeneric("modelFits"))

#' Regression coefficients from full model fit
#'
#' Access the full model fitted coefficients of an \code{\linkS4class{edgeFit}} object.
#'
#' @param object \code{\linkS4class{edgeFit}}
#' 
#' @usage betaCoef(object)
#' 
#' @return \code{betaCoef} returns the regression coefficients.
#'
#' @author John Storey, Andrew Bass 
#' 
#' @seealso \code{\link{edgeFit}}
#' 
#' @keywords betaCoef
#' 
#' @rdname betaCoef-method
#' 
#' @exportMethod betaCoef
setGeneric("betaCoef", function(object) standardGeneric("betaCoef"))

#' Statistical method used in analysis
#'
#' Access the statistic type in an \code{\linkS4class{edgeFit}} object. Can either be the Optimal Discovery Procedure (odp) or the likelihood ratio test (lrt).
#'
#' @param object \code{\linkS4class{edgeFit}}
#' 
#' @usage sType(object)
#' 
#' @return \code{sType} returns the statistic type- either "odp" or "lrt".
#'
#' @author John Storey, Andrew Bass 
#' 
#' @seealso \code{\link{edgeFit}}, \code{\linkS4class{edgeFit}} and 
#' \code{\linkS4class{edgeSet}}
#' 
#' @keywords sType
#' 
#' @rdname sType-method
#' 
#' @exportMethod sType
setGeneric("sType", function(object) standardGeneric("sType"))

#' Fitted data from the full model 
#'
#' Access the fitted data from the full model in an \code{\linkS4class{edgeFit}} object.
#'
#' @param object \code{\linkS4class{edgeFit}}
#'  
#' @usage fitFull(object)
#' 
#' @return \code{fitFull} returns a matrix of fitted values from full model.
#'
#' @author John Storey, Andrew Bass
#' 
#' @seealso \code{\link{edgeFit}}
#' 
#' @keywords fitFull
#' 
#' @rdname fitFull-method
#' 
#' @exportMethod fitFull
setGeneric("fitFull", function(object) standardGeneric("fitFull"))

#' Fitted data from the null model
#'
#' Access the fitted data from the null model in an \code{\linkS4class{edgeFit}} object.
#'
#' @param object \code{\linkS4class{edgeFit}}
#' 
#' @usage fitNull(object)
#' 
#' @return \code{fitFull} returns a matrix of fitted values from null model.
#'
#' @author  John Storey, Andrew Bass
#' 
#' @seealso \code{\link{edgeFit}}
#' 
#' @keywords fitNull
#' 
#' @rdname fitNull-method
#' 
#' @exportMethod fitNull  
setGeneric("fitNull", function(object) standardGeneric("fitNull"))

#' Residuals of full model fit
#'
#' Access the fitted full model residuals in an \code{\linkS4class{edgeFit}} object.
#'
#' @param object \code{\linkS4class{edgeFit}}
#' 
#' @usage resFull(object)
#' 
#' @return \code{resFull} returns a matrix of residuals from full model.
#'
#' @author John Storey, Andrew Bass
#' 
#' @seealso \code{\link{edgeFit}}
#' 
#' @keywords resFull
#' 
#' @rdname resFull-method
#' 
#' @exportMethod resFull
setGeneric("resFull", function(object) standardGeneric("resFull"))

#' Residuals of null model fit
#'
#' Access the fitted null model residuals in an \code{\linkS4class{edgeFit}} object.
#'
#' @param object \code{\linkS4class{edgeFit}}
#' 
#' @usage resNull(object)
#' 
#' @return \code{resNull} returns a matrix of residuals from null model.
#'
#' @author John Storey, Andrew Bass
#' 
#' @seealso  \code{\link{edgeFit}}
#' 
#' @keywords  resNull
#' 
#' @rdname resNull-method
#' 
#' @exportMethod resNull 
setGeneric("resNull", function(object) standardGeneric("resNull"))
#' Summary of edge objects
#'
#' Summary of edgeFit and edgeSet objects
#' 
#' @param object \code{\linkS4class{edgeSet}}
#' @param \dots additional parameters
#' 
#' @author John Storey, Andrew Bass
#' 
#' @keywords summary
#'
#' @export summary
setGeneric("summary")