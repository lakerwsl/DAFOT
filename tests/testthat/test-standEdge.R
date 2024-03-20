test_that("standEdge returns expected results on known input abundances", {
  P <- structure(c(0.4, 0, 0.4, 0, 0.2, 0, 0, 0, 0, 0, 0.4, 0.2, 0.4, 0, 0, 0, 0, 0, 0.8, 0, 0, 0.2, 0, 0, 0, 0, 0, 0.8, 0, 0, 0, 0.2, 0, 0, 0, 0, 0.4, 0.4, 0, 0.2, 0, 0, 0, 0, 0), dim = c(9L, 5L))
  Q <- structure(c(0.4, 0, 0.2, 0, 0.4, 0, 0, 0, 0, 0, 0.4, 0.2, 0.4, 0, 0, 0, 0, 0, 0.4, 0.2, 0.2, 0.2, 0, 0, 0, 0, 0, 0.8, 0, 0.2, 0, 0, 0, 0, 0, 0, 0.4, 0.2, 0.2, 0.2, 0, 0, 0, 0, 0), dim = c(9L, 5L))

  expect_equal(standEdge(P, Q, 5, 5), c(0.4, 0, 1.0, 0, 0, 0, 0, 0, 0), tolerance = 0.1)
})

test_that("standEdge returns no differences on the same input abundances", {
  P <- Q <- structure(c(0.4, 0, 0.4, 0, 0.2, 0, 0, 0, 0, 0, 0.4, 0.2, 0.4, 0, 0, 0, 0, 0, 0.8, 0, 0, 0.2, 0, 0, 0, 0, 0, 0.8, 0, 0, 0, 0.2, 0, 0, 0, 0, 0.4, 0.4, 0, 0.2, 0, 0, 0, 0, 0), dim = c(9L, 5L))

  expect_equal(standEdge(P, Q, 5, 5), rep(0, 9))
})
