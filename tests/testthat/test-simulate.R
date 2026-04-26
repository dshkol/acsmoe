test_that("simulation is deterministic with a fixed seed", {
  set.seed(42)
  a <- acs_simulate(c(x = 100, y = 50), c(10, 5), n_sims = 5)
  set.seed(42)
  b <- acs_simulate(c(x = 100, y = 50), c(10, 5), n_sims = 5)
  expect_equal(a, b)
})

test_that("simulation function summarizes derived values", {
  set.seed(42)
  out <- acs_simulate_fn(c(100, 50), c(10, 5), fn = sum, n_sims = 200)
  expect_named(out, c("estimate", "se"))
  expect_true(is.finite(out$estimate))
})

test_that("normal simulations recover input moments approximately", {
  set.seed(42)
  estimates <- c(100, 50)
  moes <- c(10, 5)
  draws <- acs_simulate(estimates, moes, n_sims = 50000, dist = "normal")
  expect_equal(colMeans(draws), estimates, tolerance = 0.1)
  expect_equal(apply(draws, 2, stats::sd), moe_to_se(moes), tolerance = 0.1)
})

test_that("simulate_fn forwards distribution and confidence settings", {
  set.seed(42)
  expect_warning(
    out <- acs_simulate_fn(
      estimates = c(-1, 10),
      moes = c(1, 1),
      fn = sum,
      n_sims = 100,
      dist = "truncated_normal",
      conf = 0.95
    ),
    "does not preserve"
  )
  expect_true(is.finite(out$estimate))
})
