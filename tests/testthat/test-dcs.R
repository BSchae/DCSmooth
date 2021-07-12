################################################################################
#                                                                              #
#                               Test for dcs()                                 #
#                                                                              #
################################################################################

### Tests for variance estimation
context("Variance Estimation")

test_that("model order is actually used", {
  Y = y.norm1 + rnorm(101^2)
  model_order_test = list(ar = c(2, 1), ma = c(1, 3))
  dcs_iid = dcs(Y, set.options(var_est = "iid"), model_order = model_order_test)
  dcs_qarma = dcs(Y, set.options(var_est = "qarma"), 
                  model_order = model_order_test)
  dcs_sarma = dcs(Y, set.options(var_est = "sarma"), 
                  model_order = model_order_test)
  dcs_lm = dcs(Y, set.options(var_est = "lm"), model_order = model_order_test)
  
  expect_equal(dim(dcs_iid$var_model$ar), NULL)
  expect_equal(dim(dcs_iid$var_model$ma), NULL)
  expect_equal(dim(dcs_qarma$var_model$ar), c(3, 2))
  expect_equal(dim(dcs_qarma$var_model$ma), c(2, 4))
  expect_equal(dim(dcs_sarma$var_model$ar), c(3, 2))
  expect_equal(dim(dcs_sarma$var_model$ma), c(2, 4))
  expect_equal(dim(dcs_lm$var_model$ar), c(3, 2))
  expect_equal(dim(dcs_lm$var_model$ma), c(2, 4))
})
