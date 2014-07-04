\name{gibson}
\docType{data}
\alias{gibson}

\title{Gene expression dataset from Idaghdour et al. (2008)}

\usage{
data(gibson)
}

\description{
  The data provide gene expression measurements in peripheral blood leukocyte samples from three 
  Moroccan Amazigh groups leading distinct ways of life: desert nomadic (DESERT), mountain agrarian
  (VILLAGE), and coastal urban (AGADIR).
}

\value{
A list called \code{gibson} containing:
  \item{covar}{A 4 rows by 46 columns data frame describing the study design.}
  \item{exprdat}{A 2000 rows by 46 columns matrix of gene expression values.}
}

\references{
  Idaghdour Y, Storey JD, Jadallah S, and Gibson G. (2008) A genome-wide gene 
  expression signature of lifestyle in peripheral blood of Moroccan Amazighs. 
  PLoS Genetics, 4: e1000052.
}

\examples{
data(gibson)
require(splines)
# Create an ExpressionSet
covar = data.frame(t(gibson$covar))
expSet <- ExpressionSet(assayData = gibson$exprdat, phenoData = as(covar, "AnnotatedDataFrame"))

# Create edgeSet for experiment- static experiment
mNull <- ~Gender + Batch
mFull <- ~Gender + Batch + Location
edgeObj <- edgeSet(expSet, full.model = mFull, null.model = mNull)

# Run either ODP/lrt statistic to determine significant genes
edge.odp <- odp(edgeObj, bs.its=10)
edge.lrt <- lrt(edgeObj, nullDistn="normal")

summary(edge.odp) 

}

\note{
These data are a random subset of 2000 genes from the total number of genes in the original data set.  
To download the full data set, go to \url{http://genomine.org/edge/}.
}

\keyword{datasets}  
