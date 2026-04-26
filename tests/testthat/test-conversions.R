test_that("MOE and SE conversions round-trip", {
  moe <- c(10, 25, 100)
  se <- moe_to_se(moe)
  expect_equal(se_to_moe(se), moe)
})

test_that("confidence intervals widen at higher confidence", {
  ci90 <- moe_ci(100, 10, conf_in = 0.90, conf_out = 0.90)
  ci95 <- moe_ci(100, 10, conf_in = 0.90, conf_out = 0.95)
  expect_lt(ci95$lower, ci90$lower)
  expect_gt(ci95$upper, ci90$upper)
})

test_that("CV reliability thresholds classify estimates", {
  expect_equal(as.character(acs_reliability(1000, 10)), "reliable")
  expect_equal(as.character(acs_reliability(100, 80)), "unreliable")
})

test_that("reliability passes confidence level through to CV calculation", {
  moe95 <- se_to_moe(39, conf = 0.95)
  expect_equal(as.character(acs_reliability(100, moe95, conf = 0.95)), "caveat")
  expect_equal(as.character(acs_reliability(100, moe95, conf = 0.90)), "unreliable")
})
