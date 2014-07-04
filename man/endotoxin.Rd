\name{endotoxin}
\docType{data}
\alias{endotoxin}

\title{Gene expression dataset from Calvano et al. (2005) Nature}

\usage{
  data(endotoxin)
}

\description{
  The data provide gene expression measurements in an endotoxin study where four subjects were given endotoxin and four subjects were given a placebo. Blood samples were collected and leukocytes were isolated from the samples before infusion and at times 2, 4, 6, 9, 24 hours.  
}

\value{
A list called \code{endotoxin} containing:
\item{expr}{A 5000 rows by 46 columns data frame containing expression values.}
\item{class}{A vector of length 46 containing information about which individuals
were given endotoxin}
\item{individual}{A vector of length 46 providing indexing measurements for each
individual in the experiment.}
\item{time}{A vector of length 46 indicating time measurements.}
}

\references{
    Calvano, S. E., Xiao, W., Richards, D. R., Felciano, R. M., Baker, H. V., Cho, R.J., Chen, R. O., Brownstein, B. H., Cobb, J. P., Tschoeke, S. K., et al. "A network-based analysis of systemic inflammation in humans." Nature. 2005 Oct 13; 437(7061):1032-7. Epub 2005 Aug 31. 
    
    Storey JD, Xiao W, Leek JT, Tompkins RG, and Davis RW. (2005) Significance 
    analysis of time course microarray experiments. PNAS, 102: 12837-12842.
}

\examples{
require(splines)
data(endotoxin)
# Create ExpressionSet object
expr <- endotoxin$expr
colnames(expr) <- NULL
expSet <- ExpressionSet(assayData = expr, phenoData = as(data.frame(individual=endotoxin$individual, time=endotoxin$time, class=endotoxin$class), "AnnotatedDataFrame"))

# Formulate null and full models in experiement ... Note interaction term is a way of taking into account group effects
mNull <- ~-1 + ns(time, df=4, intercept = FALSE)
mFull <- ~-1 + ns(time, df=4, intercept = FALSE) + ns(time, df=4, intercept = FALSE):class + class

# Create edgeSet object
edgeObj <- edgeSet(expSet, full.model = mFull, null.model = mNull, individual=endotoxin$individual)

# Perform ODP/lrt statistic to determine significant genes in study
edge.odp <- odp(edgeObj, bs.its=10)
edge.lrt <- lrt(edgeObj, nullDistn="bootstrap", bs.its=10)

summary(edge.odp) 
}

\note{
  The data is a random subset of 5000 genes from the full dataset. To download the full data set, go to \url{http://genomine.org/edge/}.
}

\keyword{datasets}  