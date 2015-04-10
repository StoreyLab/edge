#' Performs F-test (likelihood ratio test using Normal likelihood)
#'
#' \code{lrt} performs a generalized likelihood ratio test using the full and 
#' null models. 
#' 
#' @param object \code{S4 object}: \code{\linkS4class{deSet}}.
#' @param de.fit \code{S4 object}: \code{\linkS4class{deFit}}. Optional.
#' @param nullDistn \code{character}: either "normal" or "bootstrap", If 
#' "normal" then the p-values are calculated using the F distribution. If 
#' "bootstrap" then a bootstrap algorithm is implemented to simulate 
#' statistics from the null distribution. In the "bootstrap" case, empirical 
#' p-values are calculated using the observed and null statistics (see 
#' \code{\link{empPvals}}). Default is "normal".
#' @param bs.its \code{integer}: number of null statistics generated (only 
#' applicable for "bootstrap" method). Default is 100.
#' @param seed \code{integer}: set the seed value. Default is NULL.
#' @param verbose \code{boolean}: print iterations for bootstrap method. 
#' Default is TRUE.
#' @param ... Additional arguments for \code{\link{apply_qvalue}} and 
#' \code{\link{empPvals}} function.
#'
#' @details 
#' \code{lrt} fits the full and null models to each gene using the function
#' \code{\link{fit_models}} and then performs a likelihood ratio test. The 
#' user has the option to calculate p-values from either the F distribution or 
#' through a bootstrap algorithm. If \code{nullDistn} is "bootstrap"  
#' then empirical p-values will be determined from the \code{\link{qvalue}} 
#' package (see \code{\link{empPvals}}). 
#'
#' @author John Storey, Andrew Bass
#'
#' @return \code{\linkS4class{deSet}} object
#'
#' @examples
#' # import data
#' library(splines)
#' data(kidney)
#' age <- kidney$age
#' sex <- kidney$sex
#' kidexpr <- kidney$kidexpr
#' cov <- data.frame(sex = sex, age = age)
#' 
#' # create models 
#' null_model <- ~sex
#' full_model <- ~sex + ns(age, df = 4)
#' 
#' # create deSet object from data
#' de_obj <- build_models(data = kidexpr, cov = cov, null.model = null_model, 
#' full.model = full_model)
#'
#' # lrt method
#' de_lrt <- lrt(de_obj, nullDistn = "normal")
#'
#' # to generate p-values from bootstrap
#' de_lrt <- lrt(de_obj, nullDistn = "bootstrap", bs.its = 30)
#'
#' # input a deFit object directly
#' de_fit <- fit_models(de_obj, stat.type = "lrt")
#' de_lrt <- lrt(de_obj, de.fit = de_fit)
#' 
#' # summarize object
#' summary(de_lrt)
#'
#' @references
#' Storey JD, Xiao W, Leek JT, Tompkins RG, and Davis RW. (2005) Significance 
#' analysis of time course microarray experiments. Proceedings of the National 
#' Academy of Sciences, 102: 12837-12842.
#'
#' @seealso \code{\linkS4class{deSet}}, \code{\link{build_models}}, 
#' \code{\link{odp}}
#'
#' @export
setGeneric("lrt", function(object, de.fit,
                           nullDistn = c("normal","bootstrap"), bs.its = 100,
                           seed = NULL, verbose = TRUE, ...) 
  standardGeneric("lrt"))


#' The optimal discovery procedure  
#'
#' \code{odp} performs the optimal discovery procedure, which is a new 
#' approach for optimally performing many hypothesis tests in a 
#' high-dimensional study. When testing whether a feature is significant, the 
#' optimal discovery procedure uses information across all features when 
#' testing for significance. 
#'
#' @param object \code{S4 object}: \code{\linkS4class{deSet}}
#' @param de.fit \code{S4 object}: \code{\linkS4class{deFit}}. Optional.
#' @param odp.parms \code{list}: parameters for each cluster. See 
#' \code{\link{kl_clust}}.
#' @param bs.its \code{numeric}: number of null bootstrap iterations. Default
#' is 100.
#' @param n.mods \code{integer}: number of clusters used in 
#' \code{\link{kl_clust}}. Default is 50.
#' @param seed \code{integer}: set the seed value. Default is NULL.
#' @param verbose \code{boolean}: print iterations for bootstrap method. 
#' Default is TRUE.
#' @param ... Additional arguments for \code{\link{qvalue}} and 
#' \code{\link{empPvals}}.
#'
#'
#' @details 
#' The full ODP estimator computationally grows quadratically with 
#' respect to the number of genes. This becomes computationally infeasible at 
#' a certain point. Therefore, an alternative method called mODP is used which
#' has been shown to provide results that are very similar. mODP utilizes a 
#' k-means clustering algorithm where genes are assigned to a cluster based on 
#' the Kullback-Leiber distance. Each gene is assigned an module-average 
#' parameter to calculate the ODP score and it reduces the computations time 
#' to linear (See Woo, Leek and Storey 2010). If the number of clusters is equal 
#' to the number of genes then the original ODP is implemented. Depending on 
#' the number of hypothesis tests, this can take some time.
#'
#' @return \code{\linkS4class{deSet}} object
#'
#' @examples
#' # import data
#' library(splines)
#' data(kidney)
#' age <- kidney$age
#' sex <- kidney$sex
#' kidexpr <- kidney$kidexpr
#' cov <- data.frame(sex = sex, age = age)
#' 
#' # create models 
#' null_model <- ~sex
#' full_model <- ~sex + ns(age, df = 4)
#' 
#' # create deSet object from data
#' de_obj <- build_models(data = kidexpr, cov = cov, 
#' null.model = null_model, full.model = full_model)
#'
#' # odp method
#' de_odp <- odp(de_obj, bs.its = 30)
#'
#' # input a deFit object or ODP parameters ... not necessary
#' de_fit <- fit_models(de_obj, stat.type = "odp")
#' de_clust <- kl_clust(de_obj, n.mods = 10)
#' de_odp <- odp(de_obj, de.fit = de_fit, odp.parms = de_clust, 
#' bs.its = 30)
#' 
#' # summarize object
#' summary(de_odp)
#'
#' @references
#' Storey JD. (2007) The optimal discovery procedure: A new approach to 
#' simultaneous significance testing. Journal of the Royal Statistical 
#' Society, Series B, 69: 347-368.
#'
#' Storey JD, Dai JY, and Leek JT. (2007) The optimal discovery procedure for 
#' large-scale significance testing, with applications to comparative 
#' microarray experiments. Biostatistics, 8: 414-432.
#'
#' Woo S, Leek JT, Storey JD (2010) A computationally efficient modular 
#' optimal discovery procedure. Bioinformatics, 27(4): 509-515.
#'
#' @author John Storey, Jeffrey Leek, Andrew Bass
#'
#' @seealso \code{\link{kl_clust}}, \code{\link{build_models}} and 
#' \code{\linkS4class{deSet}}
#'
#' @export
setGeneric("odp", function(object, de.fit, odp.parms = NULL, bs.its = 100,
                           n.mods = 50, seed = NULL, verbose = TRUE, ...)
  standardGeneric("odp"))


#' Modular optimal discovery procedure (mODP) 
#'
#' \code{kl_clust} is an implementation of mODP that assigns genes to modules 
#' based off of the Kullback-Leibler distance.
#'
#' @param object \code{S4 object}: \code{\linkS4class{deSet}}.
#' @param de.fit \code{S4 object}: \code{\linkS4class{deFit}}.
#' @param n.mods \code{integer}: number of clusters.
#'
#' @details mODP utilizes a k-means clustering algorithm where genes are 
#' assigned to a cluster based on the Kullback-Leiber distance. Each gene is 
#' assigned an module-average parameter to calculate the ODP score (See Woo, 
#' Leek and Storey 2010 for more details). The mODP and full ODP produce near 
#' exact results but mODP has the advantage of being computationally
#' feasible.
#'
#' @note The results are generally insensitive to the number of modules after a
#' certain threshold of about n.mods>=50. It is recommended that users 
#' experiment with the number of clusters. If the number of clusters is equal 
#' to the number of genes then the original ODP is implemented. Depending on 
#' the number of hypothesis tests, this can take some time.
#' 
#' @return
#' A list with the following slots:
#' \itemize{
#'   \item {mu.full: cluster means from full model.}
#'   \item {mu.null: cluster means from null model.}
#'   \item {sig.full: cluster standard deviations from full model.}
#'   \item {sig.null: cluster standard deviations from null model.}
#'   \item {n.per.mod: total members in each cluster.}
#'   \item {clustMembers: cluster membership for each gene.}
#' }
#'
#' @examples
#' # import data
#' library(splines)
#' data(kidney)
#' age <- kidney$age
#' sex <- kidney$sex
#' kidexpr <- kidney$kidexpr
#' cov <- data.frame(sex = sex, age = age)
#' 
#' # create models 
#' null_model <- ~sex
#' full_model <- ~sex + ns(age, df = 4)
#' 
#' # create deSet object from data
#' de_obj <- build_models(data = kidexpr, cov = cov, null.model = null_model, 
#' full.model = full_model)
#'
#' # mODP method
#' de_clust <- kl_clust(de_obj)
#'
#' # change the number of clusters
#' de_clust <- kl_clust(de_obj, n.mods = 10)
#'
#' # input a deFit object 
#' de_fit <- fit_models(de_obj, stat.type = "odp")
#' de_clust <- kl_clust(de_obj, de.fit = de_fit)
#'
#' @references
#' Storey JD. (2007) The optimal discovery procedure: A new approach to 
#' simultaneous significance testing. Journal of the Royal Statistical 
#' Society, Series B, 69: 347-368.
#'
#' Storey JD, Dai JY, and Leek JT. (2007) The optimal discovery procedure for 
#' large-scale significance testing, with applications to comparative 
#' microarray experiments. Biostatistics, 8: 414-432.
#'
#' Woo S, Leek JT, Storey JD (2010) A computationally efficient modular optimal
#'  discovery procedure. Bioinformatics, 27(4): 509-515.
#'
#' @author John Storey, Jeffrey Leek
#'
#' @seealso \code{\link{odp}}, \code{\link{lrt}} and 
#' \code{\link{fit_models}}
#'
#' @exportMethod kl_clust
setGeneric("kl_clust", function(object, de.fit = NULL, n.mods = 50)
  standardGeneric("kl_clust"))

#' Linear regression of the null and full models 
#'
#' \code{fit_models} fits a linear model to each gene by using the least 
#' squares method. Model fits can be either statistic type "odp" (optimal 
#' discovery procedure) or "lrt" (likelihood ratio test).
#'
#' @param object \code{S4 object}: \code{\linkS4class{deSet}}.
#' @param stat.type \code{character}: type of statistic to be used. Either 
#' "lrt" or "odp". Default is "lrt".
#'
#' @details 
#' If "odp" method is implemented then the null model is removed from the full 
#' model (see Storey 2007). 
#'
#' @note \code{fit_models} does not have to be called by the user to use
#' \code{\link{odp}}, \code{\link{lrt}} or \code{\link{kl_clust}} as it is an 
#' optional input and is implemented in the methods. The 
#' \code{\linkS4class{deFit}} object can be created by the user if a different 
#' statistical implementation is required.
#'
#' @return \code{\linkS4class{deFit}} object
#'
#' @examples
#' # import data
#' library(splines)
#' data(kidney)
#' age <- kidney$age
#' sex <- kidney$sex
#' kidexpr <- kidney$kidexpr
#' cov <- data.frame(sex = sex, age = age)
#' 
#' # create models 
#' null_model <- ~sex
#' full_model <- ~sex + ns(age, df = 4)
#' 
#' # create deSet object from data
#' de_obj <- build_models(data = kidexpr, cov = cov, null.model = null_model, 
#' full.model = full_model)
#'
#' # retrieve statistics from linear regression for each gene
#' fit_lrt <- fit_models(de_obj, stat.type = "lrt") # lrt method
#' fit_odp <- fit_models(de_obj, stat.type = "odp") # odp method
#' 
#' # summarize object
#' summary(fit_odp)
#'
#' @references
#' Storey JD. (2007) The optimal discovery procedure: A new approach to 
#' simultaneous significance testing. Journal of the Royal Statistical 
#' Society, Series B, 69: 347-368.
#'
#' Storey JD, Dai JY, and Leek JT. (2007) The optimal discovery procedure for 
#' large-scale significance testing, with applications to comparative 
#' microarray experiments. Biostatistics, 8: 414-432.
#'
#' Storey JD, Xiao W, Leek JT, Tompkins RG, and Davis RW. (2005) Significance 
#' analysis of time course microarray experiments. Proceedings of the National 
#' Academy of Sciences, 102: 12837-12842.
#'
#' @seealso \code{\linkS4class{deFit}}, \code{\link{odp}} and
#' \code{\link{lrt}}
#'
#' @author John Storey
#' @exportMethod fit_models
setGeneric("fit_models",
           function(object, stat.type = c("lrt", "odp")) {
             standardGeneric("fit_models")
           })

#' Create a deSet object from an ExpressionSet
#'
#' Creates a \code{\linkS4class{deSet}} object that extends the 
#' \code{\link{ExpressionSet}} object.
#'
#' @param object \code{S4 object}: \code{\link{ExpressionSet}}
#' @param full.model \code{formula}: full model containing the both the 
#' adjustment and the biological variables for the experiment.
#' @param null.model \code{formula}: null model containing the adjustment 
#' variables for the experiment.
#' @param individual \code{factor}: information on repeated samples in 
#' experiment. 
#'
#' @note It is essential that the null and full models have the same variables
#' as the ExpressionSet phenoType column names.
#'
#' @return \code{\linkS4class{deSet}} object
#'
#' @examples
#' # import data
#' library(splines)
#' data(kidney)
#' age <- kidney$age
#' sex <- kidney$sex
#' kidexpr <- kidney$kidexpr
#' cov <- data.frame(sex = sex, age = age) 
#' pDat <- as(cov, "AnnotatedDataFrame")
#' exp_set <- ExpressionSet(assayData = kidexpr, phenoData = pDat)
#' 
#' # create models 
#' null_model <- ~sex
#' full_model <- ~sex + ns(age, df = 4)
#' 
#' # create deSet object from data
#' de_obj <- deSet(exp_set, null.model = null_model, 
#' full.model = full_model)
#'
#' # optionally add individuals to experiment, in this case there are 36 
#' # individuals that were sampled twice
#' indSamples <- as.factor(rep(1:36, each = 2))
#' de_obj <- deSet(exp_set, null.model = null_model, 
#' full.model = full_model, ind = indSamples)
#' summary(de_obj)
#' @seealso \code{\linkS4class{deSet}}, \code{\link{odp}} and 
#' \code{\link{lrt}}
#'
#' @author John Storey, Andrew Bass
#'
#' @export
setGeneric("deSet", function(object, full.model, null.model, 
                             individual=NULL) standardGeneric("deSet"))

#' Estimate the q-values for a given set of p-values
#'
#' Runs \code{\link{qvalue}} on a \code{\linkS4class{deSet}} object.
#'
#' @param object \code{S4 object}: \code{\linkS4class{deSet}}
#' @param ... Additional arguments for \code{\link{qvalue}}
#'
#' @return \code{\linkS4class{deSet}} object
#'
#' @examples
#' # import data
#' library(splines)
#' data(kidney)
#' age <- kidney$age
#' sex <- kidney$sex
#' kidexpr <- kidney$kidexpr
#' cov <- data.frame(sex = sex, age = age)
#' 
#' # create models 
#' null_model <- ~sex
#' full_model <- ~sex + ns(age, df = 4)
#' 
#' # create deSet object from data
#' de_obj <- build_models(data = kidexpr, cov = cov, null.model = null_model, 
#' full.model = full_model)
#'
#' # Run lrt (or odp) and apply_qvalue
#' de_lrt <- lrt(de_obj)
#' de_lrt <- apply_qvalue(de_lrt, fdr.level = 0.05, 
#' pi0.method = "bootstrap", adj=1.2)
#' summary(de_lrt)
#'
#' @references
#' Storey JD and Tibshirani R. (2003) Statistical significance for 
#' genome-wide studies. Proceedings of the National Academy of Sciences, 
#' 100: 9440-9445
#'
#' @seealso \code{\linkS4class{deSet}}, \code{\link{odp}} and 
#' \code{\link{lrt}}
#'
#' @author John Storey, Andrew Bass
#'
#' @export
setGeneric("apply_qvalue", function(object, ...)
  standardGeneric("apply_qvalue"))

#' Estimate surrogate variables 
#'
#' Runs \code{\link{sva}} on the null and full models in 
#' \code{\linkS4class{deSet}}. See \code{\link{sva}} for additional details.
#'
#' @param object \code{S4 object}: \code{\linkS4class{deSet}} 
#' @param ... Additional arguments for \code{\link{sva}}
#'
#' @return \code{\linkS4class{deSet}} object
#' 
#' @examples
#' # import data
#' library(splines)
#' data(kidney)
#' age <- kidney$age
#' sex <- kidney$sex
#' kidexpr <- kidney$kidexpr
#' cov <- data.frame(sex = sex, age = age)
#' 
#' # create models 
#' null_model <- ~sex
#' full_model <- ~sex + ns(age, df = 4)
#' 
#' # create deSet object from data
#' de_obj <- build_models(data = kidexpr, cov = cov, null.model = null_model, 
#' full.model = full_model)
#' 
#' # run surrogate variable analysis
#' de_sva <- apply_sva(de_obj)
#' 
#' # run odp/lrt with surrogate variables added
#' de_odp <- odp(de_sva, bs.its = 30)
#' summary(de_odp)
#' @seealso \code{\linkS4class{deSet}}, \code{\link{odp}} and 
#' \code{\link{lrt}}
#'
#' @references
#' Leek JT, Storey JD (2007) Capturing Heterogeneity in Gene Expression 
#' Studies by Surrogate Variable Analysis. PLoS Genet 3(9): e161. 
#' doi:10.1371/journal.pgen.0030161
#'
#' @author John Storey, Jeffrey Leek, Andrew Bass
#' @export
setGeneric("apply_sva", function(object, ...)
  standardGeneric("apply_sva"))

#' Supervised normalization of data in edge
#'
#' Runs \code{snm} on a deSet object based on the null and full models in 
#' \code{\linkS4class{deSet}}. See \code{\link{snm}} for additional details 
#' on the algorithm.
#'
#' @param object \code{S4 object}: \code{\linkS4class{deSet}}
#' @param int.var \code{data frame}: intensity-dependent effects.
#' @param ... Additional arguments for \code{\link{snm}}
#'
#' @return \code{apply_snm} returns an \code{\linkS4class{deSet}} object.
#'
#' @references
#' Mechan BH, Nelson PS, Storey JD. Supervised normalization of microarrays. 
#' Bioinformatics 2010;26:1308-15
#'
#' @examples
#' # simulate data
#' library(snm)
#' singleChannel <- sim.singleChannel(12345)
#' data <- singleChannel$raw.data
#' 
#' # create deSet object using build_models (can use ExpressionSet see manual)
#' cov <- data.frame(grp = singleChannel$bio.var[,2])
#' full_model <- ~grp
#' null_model <- ~1
#' 
#' # create deSet object using build_models
#' de_obj <- build_models(data = data, cov = cov, full.model = full_model, 
#' null.model = null_model)
#'
#' # run snm using intensity-dependent adjustment variable
#' de_snm <- apply_snm(de_obj, int.var = singleChannel$int.var, 
#' verbose = FALSE, num.iter = 1)
#'
#' @seealso \code{\linkS4class{deSet}}, \code{\link{odp}} and 
#' \code{\link{lrt}}
#'
#' @author John Storey, Andrew Bass
#' @export
setGeneric("apply_snm", function(object, int.var, ...) 
  standardGeneric("apply_snm"))

#' Full model equation 
#' 
#' These generic functions access and set the full model for 
#' \code{\linkS4class{deSet}} object.
#'
#' @param object \code{S4 object}: \code{\linkS4class{deSet}}
#' @param value \code{formula}: The experiment design for the full model.
#' 
#' @examples
#' # import data
#' library(splines)
#' data(kidney)
#' age <- kidney$age
#' sex <- kidney$sex
#' kidexpr <- kidney$kidexpr
#' cov <- data.frame(sex = sex, age = age)
#' 
#' # create models 
#' null_model <- ~sex
#' full_model <- ~sex + ns(age, df = 4)
#' 
#' # create deSet object from data
#' de_obj <- build_models(data = kidexpr, cov = cov, null.model = null_model, 
#' full.model = full_model)
#' 
#' # extract out the full model equation
#' mod_full <- fullModel(de_obj)
#' 
#' # change the full model in the experiment 
#' fullModel(de_obj) <- ~sex + ns(age, df = 2)
#' 
#' 
#' @return the formula for the full model.
#'    
#' @author John Storey, Andrew Bass
#' 
#' @seealso \code{\linkS4class{deSet}}
#' 
#' @export
setGeneric("fullModel", function(object) standardGeneric("fullModel"))

#' @rdname fullModel
#' @export
setGeneric("fullModel<-", function(object, value) {
  standardGeneric("fullModel<-") 
})

#' Null model equation from deSet object
#'
#' These generic functions access and set the null model for 
#' \code{\linkS4class{deSet}} object.
#'
#' @param object \code{S4 object}: \code{\linkS4class{deSet}}
#' @param value \code{formula}: The experiment design for the null model.
#' 
#' @return \code{nullModel} returns the formula for the null model.
#'
#' @examples 
#' # import data
#' library(splines)
#' data(kidney)
#' age <- kidney$age
#' sex <- kidney$sex
#' kidexpr <- kidney$kidexpr
#' cov <- data.frame(sex = sex, age = age)
#' 
#' # create models 
#' null_model <- ~sex
#' full_model <- ~sex + ns(age, df = 4)
#' 
#' # create deSet object from data
#' de_obj <- build_models(data = kidexpr, cov = cov, null.model = null_model, 
#' full.model = full_model)
#' 
#' # extract the null model equation
#' mod_null <- nullModel(de_obj)
#' 
#' # change null model in experiment but must update full model
#' nullModel(de_obj) <- ~1
#' fullModel(de_obj) <- ~1 + ns(age, df=4)
#' @author John Storey, Andrew Bass
#' 
#' @seealso \code{\linkS4class{deSet}}
#' 
#' @keywords nullModel, nullModel<-
#' 
#' @exportMethod nullModel
setGeneric("nullModel", function(object) standardGeneric("nullModel"))

#' @rdname nullModel 
#' @export
setGeneric("nullModel<-", function(object, value) {
  standardGeneric("nullModel<-") 
})

#' Matrix representation of null model
#'
#' These generic functions access and set the null matrix for 
#' \code{\linkS4class{deSet}} object.
#'
#' @param object \code{S4 object}: \code{\linkS4class{deSet}}
#' @param value \code{matrix}: null model matrix where columns are covariates 
#' and rows are observations
#' 
#' @return \code{nullMatrix} returns the value of the null model matrix.
#'
#' @examples
#' # import data
#' library(splines)
#' data(kidney)
#' age <- kidney$age
#' sex <- kidney$sex
#' kidexpr <- kidney$kidexpr
#' cov <- data.frame(sex = sex, age = age)
#' 
#' # create models 
#' null_model <- ~sex
#' full_model <- ~sex + ns(age, df = 4)
#' 
#' # create deSet object from data
#' de_obj <- build_models(data = kidexpr, cov = cov, null.model = null_model, 
#' full.model = full_model)
#' 
#' # extract the null model as a matrix
#' mat_null <- nullMatrix(de_obj)
#'
#' @author John Storey, Andrew Bass
#' 
#' @seealso \code{\linkS4class{deSet}}, \code{\link{fullModel}} and 
#' \code{\link{fullModel}}
#' 
#' @export
setGeneric("nullMatrix", function(object) standardGeneric("nullMatrix"))

#' @rdname nullMatrix
#' @export
setGeneric("nullMatrix<-", function(object, value) {
  standardGeneric("nullMatrix<-") 
})

#' Matrix representation of full model
#'
#' These generic functions access and set the full matrix for 
#' \code{\linkS4class{deSet}} object.
#'
#' @param object \code{S4 object}: \code{\linkS4class{deSet}}
#' @param value \code{matrix}: full model matrix where the columns are the 
#' covariates and rows are observations
#' 
#' @return \code{fullMatrix} returns the value of the full model matrix.
#'
#' @examples 
#' # import data
#' library(splines)
#' data(kidney)
#' age <- kidney$age
#' sex <- kidney$sex
#' kidexpr <- kidney$kidexpr
#' cov <- data.frame(sex = sex, age = age)
#' 
#' # create models 
#' null_model <- ~sex
#' full_model <- ~sex + ns(age, df = 4)
#' 
#' # create deSet object from data
#' de_obj <- build_models(data = kidexpr, cov = cov, null.model = null_model, 
#' full.model = full_model)
#' 
#' # extract the full model equation as a matrix
#' mat_full <- fullMatrix(de_obj)
#' @author Andrew Bass, John Storey
#' 
#' @seealso \code{\linkS4class{deSet}}, \code{\link{fullModel}}
#' 
#' @export
setGeneric("fullMatrix", function(object) standardGeneric("fullMatrix"))

#' @rdname fullMatrix
#' @export
setGeneric("fullMatrix<-", function(object, value) {
  standardGeneric("fullMatrix<-") 
})


#' Access/set qvalue slot
#'
#' These generic functions access and set the \code{qvalue} object in the
#' \code{\linkS4class{deSet}} object.
#' 
#' @param object \code{S4 object}: \code{\linkS4class{deSet}}
#' @param value S3 \code{object}: \code{\link{qvalue}}
#' 
#' @return  \code{qvalueObj} returns a \code{\link{qvalue}} object.
#' 
#' @examples
#' # import data
#' library(splines)
#' library(qvalue)
#' data(kidney)
#' age <- kidney$age
#' sex <- kidney$sex
#' kidexpr <- kidney$kidexpr
#' cov <- data.frame(sex = sex, age = age)
#' 
#' # create models 
#' null_model <- ~sex
#' full_model <- ~sex + ns(age, df = 4)
#' 
#' # create deSet object from data
#' de_obj <- build_models(data = kidexpr, cov = cov, null.model = null_model, 
#' full.model = full_model)
#' 
#' # run the odp method
#' de_odp <- odp(de_obj, bs.its = 20)
#' 
#' # extract out significance results
#' qval_obj <- qvalueObj(de_odp)
#' 
#' # run qvalue and assign it to deSet slot
#' pvals <- qval_obj$pvalues
#' qval_new <- qvalue(pvals, pfdr = TRUE, fdr.level = 0.1)
#' qvalueObj(de_odp) <- qval_new
#'  
#' @author John Storey, Andrew Bass
#' 
#' @seealso \code{\link{lrt}}, \code{\link{odp}} and 
#' \code{\linkS4class{deSet}}
#' 
#' @export                
setGeneric("qvalueObj", function(object) standardGeneric("qvalueObj"))

#' @rdname qvalueObj
#' @export
setGeneric("qvalueObj<-", function(object, value) {
  standardGeneric("qvalueObj<-") 
})                      

#' Individuals sampled in experiment
#'
#' These generic functions access and set the individual slot in 
#' \code{\linkS4class{deSet}}.
#' 
#' @param object \code{\linkS4class{deSet}}
#' @param value \code{factor}: identifier for each observation. Important 
#' if the same individuals are sampled multiple times.
#' 
#' @return \code{individual} returns information regarding individuals 
#' in the experiment.
#'    
#' @examples
#' library(splines)
#' # import data
#' data(endotoxin)
#' ind <- endotoxin$ind
#' time <- endotoxin$t
#' class <- endotoxin$class
#' endoexpr <- endotoxin$endoexpr
#' cov <- data.frame(individual = ind, time = time, class = class)
#' 
#' # create ExpressionSet object
#' pDat <- as(cov, "AnnotatedDataFrame")
#' exp_set <- ExpressionSet(assayData = endoexpr, phenoData = pDat)
#' 
#' # formulate null and full models in experiement
#' # note: interaction term is a way of taking into account group effects
#' mNull <- ~ns(time, df=4, intercept = FALSE)
#' mFull <- ~ns(time, df=4, intercept = FALSE) + 
#' ns(time, df=4, intercept = FALSE):class + class
#' 
#' # create deSet object
#' de_obj <- deSet(exp_set, full.model = mFull, null.model = mNull, 
#' individual = ind)
#' 
#' # extract out the individuals factor
#' ind_exp <- individual(de_obj)
#' 
#' @author John Storey, Andrew Bass
#' 
#' @seealso \code{\linkS4class{deSet}}
#' 
#' @export
setGeneric("individual", function(object) standardGeneric("individual"))

#' @rdname individual
#' @export
setGeneric("individual<-", function(object, value) {
  standardGeneric("individual<-") 
})  

#' Regression coefficients from full model fit
#'
#' Access the full model fitted coefficients of a
#' \code{\linkS4class{deFit}} object.
#'
#' @param object \code{S4 object}: \code{\linkS4class{deFit}}
#' 
#' @return \code{betaCoef} returns the regression coefficients.
#'
#' @author John Storey, Andrew Bass 
#' 
#' @seealso \code{\link{fit_models}}
#' 
#' @examples
#' # import data
#' library(splines)
#' data(kidney)
#' age <- kidney$age
#' sex <- kidney$sex
#' kidexpr <- kidney$kidexpr
#' cov <- data.frame(sex = sex, age = age)
#' 
#' # create models 
#' null_model <- ~sex
#' full_model <- ~sex + ns(age, df = 4)
#' 
#' # create deSet object from data
#' de_obj <- build_models(data = kidexpr, cov = cov, null.model = null_model, 
#' full.model = full_model)
#' 
#' # run fit_models to get model fits
#' de_fit <- fit_models(de_obj)
#' 
#' # extract beta coefficients
#' beta <- betaCoef(de_fit)
#' 
#' @export
setGeneric("betaCoef", function(object) standardGeneric("betaCoef"))

#' Statistical method used in analysis
#'
#' Access the statistic type in a \code{\linkS4class{deFit}} object. Can 
#' either be the optimal discovery procedure (odp) or the likelihood ratio 
#' test (lrt).
#'
#' @param object \code{S4 object}: \code{\linkS4class{deFit}}
#' 
#' @examples
#' # import data
#' library(splines)
#' data(kidney)
#' age <- kidney$age
#' sex <- kidney$sex
#' kidexpr <- kidney$kidexpr
#' cov <- data.frame(sex = sex, age = age)
#' 
#' # create models 
#' null_model <- ~sex
#' full_model <- ~sex + ns(age, df = 4)
#' 
#' # create deSet object from data
#' de_obj <- build_models(data = kidexpr, cov = cov, null.model = null_model, 
#' full.model = full_model)
#' 
#' # run fit_models to get model fits
#' de_fit <- fit_models(de_obj)
#' 
#' # extract the statistic type of model fits
#' stat_type <- sType(de_fit)
#' 
#' @return \code{sType} returns the statistic type- either "odp" or "lrt".
#'
#' @author John Storey, Andrew Bass 
#' 
#' @seealso \code{\link{fit_models}}, \code{\linkS4class{deFit}} and 
#' \code{\linkS4class{deSet}}
#' 
#' @keywords sType
#' 
#' @exportMethod sType
setGeneric("sType", function(object) standardGeneric("sType"))

#' Fitted data from the full model 
#'
#' Access the fitted data from the full model in a 
#' \code{\linkS4class{deFit}} object.
#'
#' @param object \code{S4 object}: \code{\linkS4class{deFit}}
#'  
#' @usage fitFull(object)
#' 
#' @return \code{fitFull} returns a matrix of fitted values from full model.
#'
#' @examples
#' # import data
#' library(splines)
#' data(kidney)
#' age <- kidney$age
#' sex <- kidney$sex
#' kidexpr <- kidney$kidexpr
#' cov <- data.frame(sex = sex, age = age)
#' 
#' # create models 
#' null_model <- ~sex
#' full_model <- ~sex + ns(age, df = 4)
#' 
#' # create deSet object from data
#' de_obj <- build_models(data = kidexpr, cov = cov, null.model = null_model, 
#' full.model = full_model)
#' 
#' # run fit_models to get model fits
#' de_fit <- fit_models(de_obj)
#' 
#' # extract fitted values for full model
#' fitted_full <- fitFull(de_fit)
#'
#' @author John Storey, Andrew Bass
#' 
#' @seealso \code{\link{fit_models}}
#' 
#' @export
setGeneric("fitFull", function(object) standardGeneric("fitFull"))

#' Fitted data from the null model
#'
#' Access the fitted data from the null model in an 
#' \code{\linkS4class{deFit}} object.
#'
#' @param object \code{S4 object}: \code{\linkS4class{deFit}}
#' 
#' @usage fitNull(object)
#' 
#' @return \code{fitNull} returns a matrix of fitted values from null model.
#'
#' @examples
#' # import data
#' library(splines)
#' data(kidney)
#' age <- kidney$age
#' sex <- kidney$sex
#' kidexpr <- kidney$kidexpr
#' cov <- data.frame(sex = sex, age = age)
#' 
#' # create models 
#' null_model <- ~sex
#' full_model <- ~sex + ns(age, df = 4)
#' 
#' # create deSet object from data
#' de_obj <- build_models(data = kidexpr, cov = cov, null.model = null_model, 
#' full.model = full_model)
#' 
#' # run fit_models to get model fits
#' de_fit <- fit_models(de_obj)
#' 
#' # extract fitted values from null model
#' fitted_null <- fitNull(de_fit) 
#'
#' @author  John Storey, Andrew Bass
#' 
#' @seealso \code{\link{fit_models}}
#' 
#' @export
setGeneric("fitNull", function(object) standardGeneric("fitNull"))

#' Residuals of full model fit
#'
#' Access the fitted full model residuals in an \code{\linkS4class{deFit}} 
#' object.
#'
#' @param object \code{S4 object}: \code{\linkS4class{deFit}}
#' 
#' @usage resFull(object)
#' 
#' @return \code{resFull} returns a matrix of residuals from full model.
#'
#' @examples
#' # import data
#' library(splines)
#' data(kidney)
#' age <- kidney$age
#' sex <- kidney$sex
#' kidexpr <- kidney$kidexpr
#' cov <- data.frame(sex = sex, age = age)
#' 
#' # create models 
#' null_model <- ~sex
#' full_model <- ~sex + ns(age, df = 4)
#' 
#' # create deSet object from data
#' de_obj <- build_models(data = kidexpr, cov = cov, null.model = null_model, 
#' full.model = full_model)
#'
#' # run fit_models to get model fits
#' de_fit <- fit_models(de_obj)
#' 
#' # extract out the full residuals from the model fit
#' res_full <- resFull(de_fit)
#'
#' @author John Storey, Andrew Bass
#' 
#' @seealso \code{\link{fit_models}}
#' 
#' @export
setGeneric("resFull", function(object) standardGeneric("resFull"))

#' Residuals of null model fit
#'
#' Access the fitted null model residuals in an \code{\linkS4class{deFit}} 
#' object.
#'
#' @param object \code{S4 object}: \code{\linkS4class{deFit}}
#' 
#' @usage resNull(object)
#' 
#' @return \code{resNull} returns a matrix of residuals from null model.
#'
#' @examples
#' # import data
#' library(splines)
#' data(kidney)
#' age <- kidney$age
#' sex <- kidney$sex
#' kidexpr <- kidney$kidexpr
#' cov <- data.frame(sex = sex, age = age)
#' 
#' # create models 
#' null_model <- ~sex
#' full_model <- ~sex + ns(age, df = 4)
#' 
#' # create deSet object from data
#' de_obj <- build_models(data = kidexpr, cov = cov, null.model = null_model, 
#' full.model = full_model)
#' 
#' # run fit_models to get model fits
#' de_fit <- fit_models(de_obj)
#' 
#' # extract out the null residuals from the model fits
#' res_null <- resNull(de_fit)
#' @author John Storey, Andrew Bass
#' 
#' @seealso  \code{\link{fit_models}}
#' 
#' @export 
setGeneric("resNull", function(object) standardGeneric("resNull"))

#' Summary of deFit and deSet
#'
#' Summary of \code{\linkS4class{deFit}} and \code{\linkS4class{deSet}} objects.
#' 
#' @param object \code{S4 object}: \code{\linkS4class{deSet}}
#' @param \dots additional parameters
#' 
#' @examples
#' # import data
#' library(splines)
#' data(kidney)
#' age <- kidney$age
#' sex <- kidney$sex
#' kidexpr <- kidney$kidexpr
#' cov <- data.frame(sex = sex, age = age)
#' 
#' # create models 
#' null_model <- ~sex
#' full_model <- ~sex + ns(age, df = 4)
#' 
#' # create deSet object from data
#' de_obj <- build_models(data = kidexpr, cov = cov, null.model = null_model, 
#' full.model = full_model)
#' 
#' # get summary
#' summary(de_obj)
#' 
#' # run odp and summarize
#' de_odp <- odp(de_obj, bs.its= 20)
#' summary(de_odp)
#' @author John Storey, Andrew Bass
#' 
#' @return
#' Summary of \code{\linkS4class{deSet}} object
#' 
#' @keywords summary
#'
#' @export summary
setGeneric("summary")

#' Show function for deFit and deSet
#'
#' Show function for \code{\linkS4class{deFit}} and \code{\linkS4class{deSet}} 
#' objects.
#' 
#' @param object \code{S4 object}: \code{\linkS4class{deSet}}
#' @param \dots additional parameters
#' 
#' @examples
#' # import data
#' library(splines)
#' data(kidney)
#' age <- kidney$age
#' sex <- kidney$sex
#' kidexpr <- kidney$kidexpr
#' cov <- data.frame(sex = sex, age = age)
#' 
#' # create models 
#' null_model <- ~sex
#' full_model <- ~sex + ns(age, df = 4)
#' 
#' # create deSet object from data
#' de_obj <- build_models(data = kidexpr, cov = cov, null.model = null_model, 
#' full.model = full_model)
#' 
#' # get summary
#' summary(de_obj)
#' 
#' # run odp and summarize
#' de_odp <- odp(de_obj, bs.its= 20)
#' de_odp
#' @author John Storey, Andrew Bass
#' 
#' @return
#' Nothing of interest
#' 
#' @export
setGeneric("show")
