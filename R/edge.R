#' @title
#' Extraction of Differential Gene Expression
#'
#' @description
#' The edge package implements methods for carrying out differential
#' expression analyses of genome-wide gene expression studies. Significance
#' testing using the optimal discovery procedure and generalized likelihood
#' ratio tests (equivalent to F-tests and t-tests) are implemented for general study
#' designs. Special functions are available to facilitate the analysis of
#' common study designs, including time course experiments. Other packages
#' such as snm, sva, and qvalue are integrated in edge to provide a wide range
#' of tools for gene expression analysis.
#'
#' @examples
#' \dontrun{
#' browseVignettes("edge")
#' }
#' @name edge
#' @author John Storey, Jeffrey Leek, Andrew Bass
#' @docType package
#' @import Biobase methods splines sva snm qvalue MASS
#' @useDynLib edge odpScoreCluster kldistance
NULL

#' @title Gene expression dataset from Calvano et al. (2005) Nature
#'
#' @description
#' The data provide gene expression measurements in an endotoxin study where
#' four subjects were given endotoxin and four subjects were given a placebo.
#' Blood samples were collected and leukocytes were isolated from the samples
#' before infusion and at times 2, 4, 6, 9, 24 hours.
#'
#' @usage data(endotoxin)
#' @format
#' \itemize{
#'   \item endoexpr: A 500 rows by 46 columns data frame containing expression
#'   values.
#'   \item class: A vector of length 46 containing information about which
#'   individuals were given endotoxin
#'   \item ind: A vector of length 46 providing indexing measurements for each
#'   individual in the experiment.
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
#' class <- endotoxin$class
#' time <- endotoxin$time
#' endoexpr <- endotoxin$endoexpr
#' cov <- data.frame(individual = ind, time = time, class = class)
#'
#' # formulate null and full models in experiement
#' # note: interaction term is a way of taking into account group effects
#' mNull <- ~ns(time, df=4, intercept = FALSE)
#' mFull <- ~ns(time, df=4, intercept = FALSE) +
#' ns(time, df=4, intercept = FALSE):class + class
#'
#' # create deSet object
#' de_obj <- build_models(endoexpr, cov = cov, full.model = mFull,
#' null.model = mNull, ind = ind)
#'
#' # Perform ODP/lrt statistic to determine significant genes in study
#' de_odp <- odp(de_obj, bs.its = 10)
#' de_lrt <- lrt(de_obj, nullDistn = "bootstrap", bs.its = 10)
#'
#' # summarize significance results
#' summary(de_odp)
#' @name endotoxin
#' @return endotoxin dataset
#' @docType data
#' @keywords datasets
NULL

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
#' Storey JD, Xiao W, Leek JT, Tompkins RG, and Davis RW. (2005) Significance
#' analysis of time course microarray experiments. PNAS, 102: 12837-12842. \cr
#' \url{http://www.pnas.org/content/100/16/9440.full}
#'
#' @examples
#' # import data
#' data(kidney)
#' sex <- kidney$sex
#' age <- kidney$age
#' kidexpr <- kidney$kidexpr
#'
#' # create model
#' de_obj <- build_study(data = kidexpr, adj.var = sex, tme = age,
#' sampling = "timecourse", basis.df = 4)
#'
#' # use the ODP/lrt method to determine significant genes
#' de_odp <- odp(de_obj, bs.its=10)
#' de_lrt <- lrt(de_obj, nullDistn = "bootstrap", bs.its = 10)
#'
#' # summarize significance results
#' summary(de_odp)
#' @name kidney
#' @return kidney dataset
#' @docType data
#' @keywords datasets
NULL

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
#' \url{http://genomine.org/de/}.
#'
#' @references
#' Idaghdour Y, Storey JD, Jadallah S, and Gibson G. (2008) A genome-wide gene
#' expression signature of lifestyle in peripheral blood of Moroccan Amazighs.
#' PLoS Genetics, 4: e1000052.
#'
#' @examples
#' # import
#' data(gibson)
#' batch <- gibson$batch
#' gender <- gibson$gender
#' location <- gibson$location
#' gibexpr <- gibson$gibexpr
#' cov <- data.frame(Batch = batch, Gender = gender,
#' Location = location)
#'
#' # create deSet for experiment- static experiment
#' mNull <- ~Gender + Batch
#' mFull <- ~Gender + Batch + Location
#'
#' # create deSet object
#' de_obj <- build_models(gibexpr, cov = cov, full.model = mFull,
#' null.model = mNull)
#'
#' # Perform ODP/lrt statistic to determine significant genes in study
#' de_odp <- odp(de_obj, bs.its = 10)
#' de_lrt <- lrt(de_obj, nullDistn = "bootstrap", bs.its = 10)
#'
#' # summarize significance results
#' summary(de_odp)
#' @name gibson
#' @return gibson dataset
#' @docType data
#' @keywords datasets
NULL
