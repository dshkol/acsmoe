test_that("zero-covariance formulas reduce to tidycensus", {
  testthat::skip_if_not_installed("tidycensus")

  estimates <- c(10, 20, 0, 0)
  moes <- c(2, 3, 4, 4)
  expect_equal(acs_sum(estimates, moes)$moe,
               tidycensus::moe_sum(moes, estimates))

  expect_equal(acs_ratio(10, 2, 50, 5)$moe,
               tidycensus::moe_ratio(10, 50, 2, 5))
  expect_equal(acs_prop(10, 2, 50, 5)$moe,
               tidycensus::moe_prop(10, 50, 2, 5))
  expect_equal(acs_product(10, 2, 50, 5)$moe,
               tidycensus::moe_product(10, 50, 2, 5))
})

test_that("proportion formula uses Census fallback when subtraction term is negative", {
  res <- acs_prop(num = 90, num_moe = 5, denom = 100, denom_moe = 20)
  fallback <- acs_ratio(num = 90, num_moe = 5, denom = 100, denom_moe = 20)
  expect_equal(res$moe, fallback$moe)
})

test_that("covariance affects differences with the correct sign", {
  no_cov <- acs_diff(100, 10, 40, 8, cov = 0)
  pos_cov <- acs_diff(100, 10, 40, 8, cov = 10)
  expect_lt(pos_cov$se, no_cov$se)
})

test_that("linear combinations use matrix covariance", {
  moes <- c(10, 20)
  se <- moe_to_se(moes)
  cov <- matrix(c(se[1]^2, 5, 5, se[2]^2), 2)
  res <- acs_linear(c(100, 50), moes, weights = c(1, -1), cov = cov)
  expect_equal(res$estimate, 50)
  expect_equal(res$se, sqrt(se[1]^2 + se[2]^2 - 10))
})

test_that("linear combinations match hand-computed zero-covariance result", {
  estimates <- c(100, 50, 25)
  moes <- c(10, 8, 6)
  weights <- c(2, -1, 0.5)
  res <- acs_linear(estimates, moes, weights)
  expected_se <- sqrt(sum((weights * moe_to_se(moes))^2))
  expect_equal(res$estimate, sum(weights * estimates))
  expect_equal(res$se, expected_se)
})

test_that("invalid covariance matrices are rejected", {
  bad <- matrix(c(1, 2, 2, 1), 2)
  expect_error(acs_linear(c(1, 1), c(1, 1), c(1, 1), cov = bad),
               "positive semidefinite")
})
