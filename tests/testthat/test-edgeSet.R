library(edge)
context("edgeSet object methods")

# create data composed of noise ----------------------------------------------- 
dat_noise <- matrix(rnorm(3000), ncol = 10)
cov <- data.frame(grp = c(rep(1, 5), rep(0, 5)))

# make edgeSet object ---------------------------------------------------------
edge_obj <- edgeModel(dat_noise, cov = cov, altMod = ~1 + grp, nullMod = ~1)

edge_obj <- lrt(edge_obj)

test_that("get methods", {
  expect_equal(fullModel(edge_obj), ~1 + grp)
  expect_equal(nullModel(edge_obj), ~1)
  expect_equal(individual(edge_obj), factor())
  
  expect_equal(class(qvalueObj(edge_obj)), "qvalue")
})

cov$new_grp <- 1:10
pData(edge_obj) <- cov
nullModel(edge_obj) <- ~1 + new_grp
fullModel(edge_obj) <- ~1 + grp + new_grp
mat_full <- model.matrix(~1 + grp + new_grp, cov)
mat_null <- model.matrix(~1 + new_grp, cov)
individual(edge_obj) <- as.factor(1:10)

test_that("set methods", {
  expect_equal(fullModel(edge_obj), ~1 + grp + new_grp)
  expect_equal(fullMatrix(edge_obj), mat_full)
  expect_error(fullModel(edge_obj) <- ~1 + DNE)
  
  expect_equal(nullModel(edge_obj), ~1 + new_grp)
  expect_equal(nullMatrix(edge_obj), mat_null) 
  expect_error(nullModel(edge_obj) <- ~1 + DNE)
  
  expect_equal(individual(edge_obj), as.factor(1:10))
  expect_error(individual(edge_obj) <- 1:10)
  
  expect_error(qvalueObj(edge_obj) <- 1:10)
})

