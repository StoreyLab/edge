#' @useDynLib edge odpScoreCluster
odp.score <- function(s.dat.cl, mu, sigma, null, m, n, cluster) {
  # Determines ODP score
  #
  # Args:
  #   s.dat.cl: Matrix of fitted data by full model
  #   mu: Vector means of clusters
  #   sigma: Vector of sd of clusters
  #   null: Boolean whether NULL model or not
  #   m: Number of genes
  #   n: Number of probes/arrays
  #   cluster: Vector of the number of members in each cluster
  #
  # Returns:
  #   scr: Vector of ODP score of each gene
  # Initilizations
  p <- length(sigma)

  # Call to C file to compute ODP score
  scr <- .C("odpScoreCluster",
            sumDat = as.double(s.dat.cl),
            mu = as.double(mu),
            sigma = as.double(sigma),
            m = as.integer(m),
            n = as.integer(n),
            p = as.integer(p),
            null = as.integer(null),
            cluster = as.integer(cluster),
            scr = double(m))$scr

  return(scr)
}

odpStat <- function(n.res, clustParms) {
  # Determines ODP statistic
  #
  # Args:
  #   n.res: null residuals
  #   clustParms: clustering parameters
  #
  # Returns:
  #   matrix of null statistics
  # Probabilities of alt and null distributions
  s.dat1 = c(t(n.res), t(clustParms$mu.full))
  s.dat0 = c(t(n.res), t(clustParms$mu.null))
  cl.den <- odp.score(s.dat0,
                      mu = rep(0, length(clustParms$sig.null)),
                      sigma = clustParms$sig.null,
                      null = TRUE,
                      m = nrow(n.res),
                      n = ncol(n.res),
                      cluster = clustParms$n.per.mod)
  cl.num <- odp.score(s.dat1,
                      mu = rowSums(clustParms$mu.full ^ 2),
                      sigma = clustParms$sig.full,
                      null = FALSE,
                      m = nrow(n.res),
                      n = ncol(n.res),
                      cluster = clustParms$n.per.mod)

  # ODP statistic
  odp.stat <-  2 * cl.num / (cl.den + cl.num)
  return(odp.stat)
}
