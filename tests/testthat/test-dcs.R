################################################################################
#                                                                              #
#                               Test for dcs()                                 #
#                                                                              #
################################################################################

### Test for incorrect inputs in dcs()
context("DCS input")

test_that("Y has exception handling", {
  expect_error(dcs(1), "Y must be a numeric matrix.")
  Y = matrix(rnorm(4^2), 4, 4)
  expect_error(dcs(Y), "Y has to be at least of dimension 5 in each direction.")
  Y = matrix(rnorm(100), 10, 10)
  Y[2, 5] = "test"
  expect_error(dcs(Y), "Y must be a numeric matrix.")
  Y[2, 5] = NA
  expect_error(dcs(Y), "Y contains missing values")
})

test_that("dcs_options has exception handling", {
  Y = matrix(rnorm(100), 10, 10)
  expect_error(dcs(Y, dcs_options = "test"),
               "Incorrect options specified, please use")
  dcs_options = set.options()
  dcs_options$type = "test"
  expect_error(dcs(Y, dcs_options), "Unsupported regression type.")
  dcs_options = set.options()
  dcs_options$test = 1
  expect_warning(dcs(Y, dcs_options), "unknown and will be ignored.")
  dcs_options = set.options()
  dcs_options$type = NULL
  expect_error(dcs(Y, dcs_options), "not specified.")
})

test_that("dcs_options default values are correct", {
  Y = matrix(rnorm(100), 10, 10)
  dcs_options = set.options()
  expect_equal(dcs(Y)$dcs_options, dcs_options)
})

test_that("dcs_options is correctly used", {
  Y = matrix(rnorm(100), 10, 10)
  dcs_options = set.options()
  expect_equal(dcs(Y, dcs_options)$dcs_options, dcs_options)
  dcs_options = set.options(type = "KR")
  expect_equal(dcs(Y, dcs_options)$dcs_options, dcs_options)
})

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
