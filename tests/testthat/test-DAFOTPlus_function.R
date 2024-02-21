test_that("IndDAFOT function works as expected", {
  d <- test_IndDAFOT
  res <- IndDAFOT(d$P, d$Y, d$tree)
  
  expect_true(all(res$P < 0.05))
})

test_that("ConIndDAFOT function works as expected", {
  
})
