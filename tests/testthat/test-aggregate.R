test_that("aggregation sums paired estimate and MOE columns", {
  dat <- data.frame(
    region = c("a", "a", "b"),
    total = c(10, 20, 30),
    total_moe = c(2, 3, 4)
  )

  out <- acs_aggregate(dat, "region", "total", "total_moe")
  expect_equal(out$total, c(30, 30))
  expect_equal(out$total_moe[1], acs_sum(c(10, 20), c(2, 3))$moe)
})

test_that("constant aggregation covariance is interpreted as correlation", {
  dat <- data.frame(
    region = c("a", "a"),
    total = c(10, 20),
    total_moe = c(2, 3)
  )
  zero <- acs_aggregate(dat, "region", "total", "total_moe",
                        cov_strategy = "zero")
  corr <- acs_aggregate(dat, "region", "total", "total_moe",
                        cov_strategy = "constant", cov_value = 0.5)
  expect_gt(corr$total_moe, zero$total_moe)
})
