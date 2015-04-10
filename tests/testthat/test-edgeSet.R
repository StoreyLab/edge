library(edge)
context("deSet object methods")

# create data composed of noise -----------------------------------------------
dat_noise <- matrix(rnorm(3000), ncol = 10)
cov <- data.frame(grp = c(rep(1, 5), rep(0, 5)))

# make deSet object ---------------------------------------------------------
de_obj <- build_models(dat_noise, cov = cov, full.model = ~1 + grp,
                      null.model = ~1)

de_obj <- lrt(de_obj)

test_that("get methods", {
  expect_equal(fullModel(de_obj), ~1 + grp)
  expect_equal(nullModel(de_obj), ~1)
  expect_equal(individual(de_obj), factor())

  expect_equal(class(qvalueObj(de_obj)), "qvalue")
})

cov$new_grp <- 1:10
pData(de_obj) <- cov
nullModel(de_obj) <- ~1 + new_grp
fullModel(de_obj) <- ~1 + grp + new_grp
mat_full <- model.matrix(~1 + grp + new_grp, cov)
mat_null <- model.matrix(~1 + new_grp, cov)
individual(de_obj) <- as.factor(1:10)

test_that("set methods", {
  expect_equal(fullModel(de_obj), ~1 + grp + new_grp)
  expect_equal(fullMatrix(de_obj), mat_full)
  expect_error(fullModel(de_obj) <- ~1 + DNE)

  expect_equal(nullModel(de_obj), ~1 + new_grp)
  expect_equal(nullMatrix(de_obj), mat_null)
  expect_error(nullModel(de_obj) <- ~1 + DNE)

  expect_equal(individual(de_obj), as.factor(1:10))
  expect_error(individual(de_obj) <- 1:10)

  expect_error(qvalueObj(de_obj) <- 1:10)
})

