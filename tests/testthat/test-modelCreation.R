library(edge)
library(splines)
context("Model creation: edgeSet/edgeModel/edgeStudy")

ngenes <- 100
nobs <- 20
# create data composed of noise -----------------------------------------------
dat_noise <- matrix(rnorm(ngenes*nobs), ncol = nobs)
cov <- data.frame(grp = c(rep(1, nobs/2), rep(0, nobs/2)))

# edgeModel -------------------------------------------------------------------
edge_objM <- edgeModel(dat_noise, cov = cov, full.model = ~1 + grp, null.model = ~1)
edge_objMi <- edgeModel(dat_noise, cov = cov, full.model = ~1 + grp, 
                       null.model = ~1, ind = factor(1:20))

# edgeStudy -------------------------------------------------------------------
edge_objS <- edgeStudy(dat_noise, grp = as.factor(cov$grp), sampling = "static")
edge_objSi <- edgeStudy(dat_noise, grp = as.factor(cov$grp), sampling = "static", 
                       ind = factor(1:20))
adj <- rnorm(20)
tme <- rnorm(20)
cov$adj <- adj
cov$tme <- tme
edge_objSit <- edgeStudy(dat_noise, grp = as.factor(cov$grp), 
                        adj.var = adj, sampling = "timecourse", 
                        ind = factor(1:20), tme = tme)

# edgeSet ---------------------------------------------------------------------
exp_set <- ExpressionSet(assayData = dat_noise, 
                         phenoData = as(cov, "AnnotatedDataFrame"))
edge_objE <- edgeSet(exp_set, full.model = ~1 + grp, null.model = ~1)
edge_objEi <- edgeSet(exp_set, full.model = ~1 + grp, null.model = ~1, 
                      ind = factor(1:20))


test_that("edgeModel method", {
  expect_equal(fullModel(edge_objM), ~1 + grp)
  expect_equal(nullModel(edge_objM), ~1) 
  
  expect_equivalent(fullMatrix(edge_objM), model.matrix(~1 + grp, cov))
  
  expect_equivalent(exprs(edge_objM), dat_noise)
  
  expect_equivalent(individual(edge_objMi), factor(1:20))
  
})

test_that("edgeStudy method", {
  expect_equal(fullModel(edge_objS), ~grp)
  expect_equal(nullModel(edge_objS), ~1)
  
  expect_equivalent(fullMatrix(edge_objS), model.matrix(~1 + grp, cov))
  
  expect_equivalent(exprs(edge_objS), dat_noise)
  
  expect_equivalent(individual(edge_objSi), factor(1:20))
  
  expect_equivalent(fullMatrix(edge_objSit), model.matrix(~adj + grp + ns(tme, df=2, intercept=FALSE) + grp:ns(tme, df = 2, intercept=FALSE),  cov))
  
})

test_that("edgeSet method", {
  expect_equal(fullModel(edge_objE), ~1 + grp)
  expect_equal(nullModel(edge_objE), ~1)
  
  expect_equivalent(fullMatrix(edge_objE), model.matrix(~1 + grp, cov))
  
  expect_equivalent(exprs(edge_objE), dat_noise)
  
  expect_equivalent(individual(edge_objEi), factor(1:20))
  
})