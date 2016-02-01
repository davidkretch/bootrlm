context("bootrlm")

test_that("un-bootstrapped results match MASS::rlm", {
  set.seed(1)
  res <- bootrlm(stack.loss ~ ., stackloss, r = 0, method = "MM")
  set.seed(1)
  cmp <- MASS::rlm(stack.loss ~ ., stackloss, method = "MM")

  expect_equal(as.numeric(res$coefficients),
               as.numeric(cmp$coefficients),
               tolerance = 1e-3)

  expect_equal(as.numeric(res$scale),
               as.numeric(cmp$s),
               tolerance = 1e-3)
})
