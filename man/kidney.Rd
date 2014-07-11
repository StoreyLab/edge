\name{kidney}
\docType{data}
\alias{kidney}

\title{Gene expression dataset from Rodwell et al. (2004)}

\description{
  Gene expression measurements from kidney samples were obtained from 72 human subjects ranging in age from 27 to 92 years. Only one array was obtained per sample, and the age and tissue type of each subject was recorded.
}

\usage{
data(kidney)
}

\value{
A list called \code{kidney} containing:
\item{kidcov}{A 133 rows by 6 columns data frame detailing the study design.}
\item{kidexpr}{A 1500 rows by 133 columns matrix of gene expression values, where each row corresponds to a different probe-set and each column to a different tissue sample.}
\item{age}{A vector of length 133 giving the age of each sample.}
\item{sex}{A vector of length 133 giving the sex of each sample.}
\item{tissue}{A vector of length 133 giving the tissue type of each sample.}
}

\references{
    Rodwell GE et al. (2004) A transcriptional profile of aging in the human kidney. 
    PLoS Biology, 2(12): e427.

    Storey JD, Xiao W, Leek JT, Tompkins RG, and Davis RW. (2005) Significance 
    analysis of time course microarray experiments. PNAS, 102: 12837-12842.
}


\examples{
data(kidney)

# Interested in cortex samples 
sex <- kidney$sex
age <- kidney$age
kidexpr <- kidney$kidexpr

# Create model
edge.obj <- edgeModel(dat=kidexpr, adj.var=data.frame(sex), tme=age, sampling="timecourse", basis.df=4)

# Use the ODP/lrt method to determine significant genes
edge.odp <- odp(edge.obj, bs.its=10)
edge.lrt <- lrt(edge.obj, nullDistn="bootstrap", bs.its=10)

summary(edge.lrt)
}

\note{
These data are a random subset of 1500 probe-sets from the total number of probe-sets in the original data set.  
To download the full data set, go to \url{http://genomine.org/edge/}.

The \code{age}, \code{sex}, and \code{tissue} data are contained in \code{kidcov} data frame.
}

\keyword{datasets}  
