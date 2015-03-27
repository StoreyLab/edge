#' Extraction and Analysis of Differential Gene Expression
#'
#' @name edge
#' @docType package
#' @import Biobase methods splines sva snm qvalue MASS
#' @useDynLib edge odpScoreCluster kldistance
NULL

#' @name endotoxin
#' @title Gene expression dataset from Calvano et al. (2005) Nature
#' 
#' @usage
#' data(endotoxin)
#' 
#' @description
#' The data provide gene expression measurements in an endotoxin study where 
#' four subjects were given endotoxin and four subjects were given a placebo. 
#' Blood samples were collected and leukocytes were isolated from the samples 
#' before infusion and at times 2, 4, 6, 9, 24 hours.  
#' 
#' @format
#' A list called \code{endotoxin} containing:
#' \itemize{
#'   \item endoexpr: A 500 rows by 46 columns data frame containing expression 
#'   values.
#'   \item class: A vector of length 46 containing information about which 
#'   individuals were given endotoxin
#'   \item ind: A vector of length 46 providing indexing measurements for each 
#'   \item individual: in the experiment.
#'   \item time: A vector of length 46 indicating time measurements.
#' }
#' 
#' @note
#' The data is a random subset of 500 genes from the full dataset. To 
#' download the full data set, go to \url{http://genomine.org/edge/}.
#' 
#' @references
#' Calvano, S. E., Xiao, W., Richards, D. R., Felciano, R. M., Baker, H. 
#' V., Cho, R.J., Chen, R. O., Brownstein, B. H., Cobb, J. P., Tschoeke, 
#' S. K., et al. "A network-based analysis of systemic inflammation in humans."
#'  Nature. 2005 Oct 13; 437(7061):1032-7. Epub 2005 Aug 31. 
#'  
#' Storey JD, Xiao W, Leek JT, Tompkins RG, and Davis RW. (2005) Significance 
#' analysis of time course microarray experiments. PNAS, 102: 12837-12842. \cr
#' \url{http://www.pnas.org/content/100/16/9440.full}
#' 
#' @examples
#' library(splines)
#' # import data
#' data(endotoxin)
#' ind <- endotoxin$ind
#' time <- endotoxin$time
#' endoexpr <- endotoxin$endoexpr
#' class <- endotoxin$class
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
#' # create edgeSet object
#' edge_obj <- edgeSet(exp_set, full.model = mFull, null.model = mNull, 
#' individual = ind)
#' 
#' # Perform ODP/lrt statistic to determine significant genes in study
#' edge_odp <- odp(edge_obj, bs.its = 10)
#' edge_lrt <- lrt(edge_obj, nullDistn = "bootstrap", bs.its = 10)
#' 
#' # summarize significance results
#' summary(edge_odp)
#' @docType data
#' @keywords datasets
NULL

#' @name kidney
#' @title Gene expression dataset from Rodwell et al. (2004)
#' 
#' @usage
#' data(kidney)
#' 
#' @description
#' Gene expression measurements from kidney samples were obtained from 72 
#' human subjects ranging in age from 27 to 92 years. Only one array was 
#' obtained per sample, and the age and tissue type of each subject was 
#' recorded.
#' 
#' @format 
#' A list called \code{kidney} containing:
#' \itemize{
#'   \item kidcov: A 133 rows by 6 columns data frame detailing the study 
#'   design.
#'   \item kidexpr: A 500 rows by 133 columns matrix of gene expression values, 
#'   where each row corresponds to a different probe-set and each column to a 
#'   different tissue sample.
#'   \item age: A vector of length 133 giving the age of each sample.
#'   \item sex: A vector of length 133 giving the sex of each sample.
#'   \item tissue: A vector of length 133 giving the tissue type of each 
#'   sample.
#' }
#' @note
#' These data are a random subset of 500 probe-sets from the total number of 
#' probe-sets in the original data set. To download the full data set, go to 
#' \url{http://genomine.org/edge/}. The \code{age}, \code{sex}, and
#' \code{tissue} data are contained in \code{kidcov} data frame.
#' 
#' @references
#' Rodwell GE et al. (2004) A transcriptional profile of aging in the human 
#' kidney. PLoS Biology, 2(12): e427.
#' 
#' Storey JD, Xiao W, Leek JT, Tompkins RG, and Davis RW. (2005) Significance 
#' analysis of time course microarray experiments. PNAS, 102: 12837-12842. \cr
#' \url{http://www.pnas.org/content/100/16/9440.full}
#' 
#' @examples
#' # import data
#' data(kidney)
#' 
#' #interested in cortex samples 
#' sex <- kidney$sex
#' age <- kidney$age
#' kidexpr <- kidney$kidexpr
#' 
#' # create model
#' edge_obj <- edgeStudy(data = kidexpr, adj.var = sex, tme = age, 
#' sampling = "timecourse", basis.df = 4)
#' 
#' # use the ODP/lrt method to determine significant genes
#' edge_odp <- odp(edge_obj, bs.its=10)
#' edge_lrt <- lrt(edge_obj, nullDistn = "bootstrap", bs.its = 10)
#' 
#' # summarize significance results
#' summary(edge_odp)
#' 
#' @docType data
#' @keywords datasets
NULL

#' @name gibson
#' @title Gene expression dataset from Idaghdour et al. (2008)
#' 
#' @usage
#' data(gibson)
#' 
#' @description
#' The data provide gene expression measurements in peripheral blood leukocyte 
#' samples from three Moroccan Amazigh groups leading distinct ways of life: 
#' desert nomadic (DESERT), mountain agrarian (VILLAGE), and coastal urban 
#' (AGADIR).
#' 
#' @format
#' A list called \code{gibson} containing:
#' \itemize{
#'   \item batch: Batches in experiment.
#'   \item location: Location of Moroccan Amazigh groups.
#'   \item gender: Gender of individuals.
#'   \item gibexpr: A 500 rows by 46 columns matrix of gene expression values.
#' }
#' 
#' @note
#' These data are a random subset of 500 genes from the total number of genes 
#' in the original data set. To download the full data set, go to 
#' \url{http://genomine.org/edge/}.
#' 
#' @references
#' Idaghdour Y, Storey JD, Jadallah S, and Gibson G. (2008) A genome-wide gene 
#' expression signature of lifestyle in peripheral blood of Moroccan Amazighs. 
#' PLoS Genetics, 4: e1000052.
#' 
#' @examples
#' # import
#' data(gibson)
#' 
#' # create an ExpressionSet
#' covar <- data.frame(Batch = gibson$batch, Gender = gibson$gender, 
#' Location = gibson$location)
#' pDat <- as(covar, "AnnotatedDataFrame")
#' expSet <- ExpressionSet(assayData = gibson$gibexpr, phenoData = pDat)
#' 
#' # create edgeSet for experiment- static experiment
#' mNull <- ~Gender + Batch
#' mFull <- ~Gender + Batch + Location
#' 
#' # create edgeSet object
#' edge_obj <- edgeSet(expSet, full.model = mFull, null.model = mNull)
#' 
#' # Perform ODP/lrt statistic to determine significant genes in study
#' edge_odp <- odp(edge_obj, bs.its = 10)
#' edge_lrt <- lrt(edge_obj, nullDistn = "bootstrap", bs.its = 10)
#' 
#' # summarize significance results
#' summary(edge_odp) 
#' 
#' @docType data
#' @keywords datasets
NULL