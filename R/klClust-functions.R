klmod <- function(obj.edgeFit, n.mods = 50) {
  # Clustering algorithm
  #
  # Args:
  #   fitFull: fitted data for full model
  #   resFull: residual data for full.model
  #   n.mods: number of clusters
  #   df: degree of freedom of null and full models
  #
  # Returns:
  #   vector of cluster assignments for each gene
  
  # Initializations
  m <- nrow(obj.edgeFit@fit.full)
  n <- ncol(obj.edgeFit@fit.full)
  if (m <= n.mods) {
    mod.member <- as.factor(1:m)
    return(mod.member)
  }
  sigma2 <- rowSums(obj.edgeFit@res.full ^ 2) / (n - ncol(obj.edgeFit@beta.coef))
  int.n.mods <- n.mods
  orig.n.mods <- n.mods
  
  int.center <- sample(x = m, 
                       size = n.mods,
                       replace = FALSE)    
  center.fitFull <- obj.edgeFit@fit.full[int.center, ]
  center.var <- sigma2[int.center]
  
  eps <- 0.1
  mod.member <- NULL
  KL <- matrix(nrow = m, 
               ncol = n.mods)
  itr <- 0
  KL.cutoff <- 1
  
  pos.center.fitFull <- center.fitFull
  pos.center.var <- center.var
  while (KL.cutoff > eps) {
    itr <- itr + 1
    pre.center.fitFull <- pos.center.fitFull
    pre.center.var <- pos.center.var
    
    temp.center.fitFull <- as.vector(t(center.fitFull))
    temp.fitFull <- as.vector(t(obj.edgeFit@fit.full))
    
    kldd <- t(matrix(kl(temp.center.fitFull, temp.fitFull, center.var, 
                        sigma2, n=n), ncol=m))
    mod.member = apply(kldd, 1, function(x) which.min(x))  	
    
    # First of all, we check whether there is any cluster that does not 
    # include any gene. For this case, we exclude this cluster from the 
    # original clusters. Therefore, it reduces the number of clusters 
    notempty <- 1:n.mods %in% unique(mod.member)
#    notempty <- sort(unique(mod.member))
   # all.equal(notempty, notempty2)
    center.fitFull <- center.fitFull[notempty, ]
    center.var <- center.var[notempty]
    KL <- KL[notempty, ]
    
    # Once the number of clusters were decided, we need to find new centers 
    # for each cluster 
    if (any(!notempty)) {
      n.mods <- sum(!notempty)
    }
    
    # Average the mean and variance over genes included in each cluster 
    l <- 1
    for (i in 1:orig.n.mods) { 
      ntmp <- sum(mod.member == i)
      if (ntmp == 0) {
        next
      } else {
        if (ntmp == 1) {
          center.fitFull[l, ] <- obj.edgeFit@fit.full[mod.member == i, ]
        } else {
          center.fitFull[l, ] <- colMeans(obj.edgeFit@fit.full[mod.member == i, ])
        }
        center.var[l] <- drop(sum(sigma2[mod.member == i]) / ntmp)
        l <- l + 1
      }
    }
    
    pos.center.fitFull <- center.fitFull
    pos.center.var <- center.var
    if (length(pos.center.var) != length(pre.center.var)) {
      # if the n.mods is reduced 
      KL.cutoff <- 1
    } else {
      KL.cutoff <- NULL
      res2 <- rowSums((pos.center.fitFull - pre.center.fitFull) ^ 2)
      normconst <- 1 / pos.center.var + 1 / pre.center.var
      centerconst <- n * ((pos.center.var / pre.center.var + pre.center.var / pos.center.var) / 2 - 1)
      KL.cutoff <- res2 / normconst + centerconst      
    }
    KL.cutoff <- max(KL.cutoff)
  }
  return(as.factor(mod.member))
}

# Clustering parameters
mod.parms <- function(obj.edgeFit, clMembers) {
            # Calculates parameters of interest from clusters
            #
            # Args:
            #   object: edge object (S4 class)
            #
            # Returns:
            #   list:
            #     mu.full: mean of clusters from full model
            #     mu.null: mean of clusters from null model
            #     sig.full: sd of clusters from full model
            #     sig.null: sd of clusters from null model
            #     n.per.mod: total members in each cluster
            #     clustMembers: members for each gene
            # Initlizations
            n <- ncol(obj.edgeFit@res.full)
            nf <- ncol(obj.edgeFit@beta.coef)
            varFull <- rowSums(obj.edgeFit@res.full ^ 2) / (n - nf)
            varNull <- rowSums(obj.edgeFit@res.null ^ 2) / (n - nf)
            mod.membership <- clMembers
            n.mods <- length(unique(mod.membership))
            
            mod.fitFull <- matrix(nrow = n.mods,
                                  ncol = n)
            n.per.mod <- vector(length = n.mods)
            mod.varNull <- vector(length = n.mods)
            mod.varFull <- vector(length = n.mods)
            # Calculate statistics (variance and mean) for each cluster
            for (i in 1:n.mods) {
              if(length(mod.membership[mod.membership == i]) == 1) {
                n.per.mod[i] <- 1
                mod.fitFull[i, ] <- obj.edgeFit@fit.full[mod.membership==i, ]
              } else {
                n.per.mod[i] <- sum(mod.membership == i)
                mod.fitFull[i, ] <- colMeans(obj.edgeFit@fit.full[mod.membership == i, ])
              }			
              mod.varNull[i] <- mean(varNull[mod.membership == i])
              mod.varFull[i] <- mean(varFull[mod.membership == i])
            }
            mod.fitNull <- 0*mod.fitFull
	         # Assign slots
            return(list(mu.full = mod.fitFull, sig.full = sqrt(mod.varFull), 
                        mu.null = mod.fitNull, sig.null = sqrt(mod.varNull), 
                        n.per.mod = n.per.mod, clustMembers = clMembers))
          }

kl <- function(temp.center.fitFull, temp.fitFull, center.var, sigma2, n) {
  # Calculates kl distance
  #
  # Args:
  #   temp.center.fitFull: Matrix of center locations for each cluster
  #   temp.fitFull: Matrix of fitted data from full model
  #   center.var: Vector of variances for cluster
  #   sigma2: Vector of variances for each gene
  #   n: Number of genes
  #
  # Returns:
  #   kldd: A matrix of kl distances to each cluster for each gene
  # Initializations
  m <- length(sigma2)
  n.cluster <- length(center.var)
  if (!is.loaded("kldistance")) {
    dyn.load(paste("kldistance", .Platform$dynlib.ext, sep=""))
  }
  # C function to calculate kl distance
  kldd <- .C("kldistance", 
             centerFit=as.double(temp.center.fitFull),
             centerVar=as.double(center.var),
             fit=as.double(temp.fitFull),
             var=as.double(sigma2),
             m=as.integer(m),
             nc=as.integer(n.cluster),
             n=as.integer(n),
             kldd=double(m * n.cluster))$kldd
  return(kldd)
}  
