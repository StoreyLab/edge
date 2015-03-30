# test whether lrt and odp generate uniform p-values --------------------------
library(edge)
context("Double ks-test: Uniform p-values")

pvals <- matrix(nrow = 5, ncol = 2)
nexp <- 5
ngenes <- 3000
nobs <- 24
for (i in 1:nexp) {
  set.seed(2000+i)
  # create data composed of noise ---------------------------------------------
  dat_noise <- matrix(rnorm(ngenes*nobs), ncol = nobs)
  cov <- data.frame(grp = c(rep(1, nobs/2), rep(0, nobs/2)))
  
  # make edgeSet object -------------------------------------------------------
  edge_obj <- edgeModel(dat_noise, cov = cov, altMod = ~1 + grp, nullMod = ~1)
  
  # run statistical tests and store p-values ----------------------------------
  edge_lrt <- lrt(edge_obj)
  edge_odp <- odp(edge_obj, bs.its = 100, n.mods=50)
  pvals[i,] <- cbind(ks.test(edge_lrt@qvalueObj$pvalues, "punif")$p.value,
                     ks.test(edge_odp@qvalueObj$pvalues, "punif")$p.value)
}

lrt_dks <- ks.test(pvals[, 1], "punif")$p
odp_dks <- ks.test(pvals[, 2], "punif")$p

test_that("likelihood ratio test double ks-test", {
  expect_more_than(lrt_dks, 0.05)
})

test_that("optimal discovery procedure double ks-test", {
  expect_more_than(odp_dks, 0.05)
})
