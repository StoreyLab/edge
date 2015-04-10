library(edge)
library(splines)
context("Model creation: deSet/build_models/build_study")

ngenes <- 100
nobs <- 20
# create data composed of noise -----------------------------------------------
dat_noise <- matrix(rnorm(ngenes*nobs), ncol = nobs)
cov <- data.frame(grp = c(rep(1, nobs/2), rep(0, nobs/2)))

# edgeModel -------------------------------------------------------------------
de_objM <- build_models(dat_noise, cov = cov, full.model = ~1 + grp, null.model = ~1)
de_objMi <- build_models(dat_noise, cov = cov, full.model = ~1 + grp,
                       null.model = ~1, ind = factor(1:20))

# edgeStudy -------------------------------------------------------------------
de_objS <- build_study(dat_noise, grp = as.factor(cov$grp), sampling = "static")
de_objSi <- build_study(dat_noise, grp = as.factor(cov$grp), sampling = "static",
                         ind = factor(1:20))
adj <- rnorm(20)
tme <- rnorm(20)
cov$adj <- adj
cov$tme <- tme
de_objSit <- build_study(dat_noise, grp = as.factor(cov$grp),
                        adj.var = adj, sampling = "timecourse",
                        ind = factor(1:20), tme = tme)

# deSet ---------------------------------------------------------------------
exp_set <- ExpressionSet(assayData = dat_noise,
                         phenoData = as(cov, "AnnotatedDataFrame"))
de_objE <- deSet(exp_set, full.model = ~1 + grp, null.model = ~1)
de_objEi <- deSet(exp_set, full.model = ~1 + grp, null.model = ~1,
                      ind = factor(1:20))


test_that("build_models method", {
  expect_equal(fullModel(de_objM), ~1 + grp)
  expect_equal(nullModel(de_objM), ~1)

  expect_equivalent(fullMatrix(de_objM), model.matrix(~1 + grp, cov))

  expect_equivalent(exprs(de_objM), dat_noise)

  expect_equivalent(individual(de_objMi), factor(1:20))

})

test_that("build_study method", {
  expect_equal(fullModel(de_objS), ~grp)
  expect_equal(nullModel(de_objS), ~1)

  expect_equivalent(fullMatrix(de_objS), model.matrix(~1 + grp, cov))

  expect_equivalent(exprs(de_objS), dat_noise)

  expect_equivalent(individual(de_objSi), factor(1:20))

  expect_equivalent(fullMatrix(de_objSit), model.matrix(~adj + grp + ns(tme, df=2, intercept=FALSE) + grp:ns(tme, df = 2, intercept=FALSE),  cov))

})

test_that("deSet method", {
  expect_equal(fullModel(de_objE), ~1 + grp)
  expect_equal(nullModel(de_objE), ~1)

  expect_equivalent(fullMatrix(de_objE), model.matrix(~1 + grp, cov))

  expect_equivalent(exprs(de_objE), dat_noise)

  expect_equivalent(individual(de_objEi), factor(1:20))

})
