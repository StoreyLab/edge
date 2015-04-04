# test whether lrt and odp generate uniform p-values --------------------------
library(edge)
context("Double ks-test: Uniform p-values")

pvals <- matrix(nrow = 5, ncol = 2)
nexp <- 1
ngenes <- 3000
nobs <- 24
for (i in 1:nexp) {
  set.seed(2000+i)
  # create data composed of noise ---------------------------------------------
  dat_noise <- matrix(rnorm(ngenes*nobs), ncol = nobs)
  cov <- data.frame(grp = c(rep(1, nobs/2), rep(0, nobs/2)))
  
  # make deSet object -------------------------------------------------------
  de_obj <- build_models(dat_noise, cov = cov, full.model = ~1 + grp,
                          null.model = ~1)
  
  # run statistical tests and store p-values ----------------------------------
  de_lrt <- lrt(de_obj)
  de_odp <- odp(de_obj, bs.its = 100, n.mods=50)
  pvals[i,] <- cbind(ks.test(de_lrt@qvalueObj$pvalues, "punif")$p.value,
                     ks.test(de_odp@qvalueObj$pvalues, "punif")$p.value)
}

lrt_dks <- ks.test(pvals[, 1], "punif")$p
odp_dks <- ks.test(pvals[, 2], "punif")$p

test_that("likelihood ratio test double ks-test", {
  expect_more_than(lrt_dks, 0.05)
})

test_that("optimal discovery procedure double ks-test", {
  expect_more_than(odp_dks, 0.05)
})
