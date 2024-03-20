P <- structure(c(0.4, 0, 0.4, 0, 0.2, 0, 0, 0, 0, 0, 0.4, 0.2, 0.4, 0, 0, 0, 0, 0, 0.8, 0, 0, 0.2, 0, 0, 0, 0, 0, 0.8, 0, 0, 0, 0.2, 0, 0, 0, 0, 0.4, 0.4, 0, 0.2, 0, 0, 0, 0, 0), dim = c(9L, 5L))
Y <- c(0.4, 0, 0.4, 0, 0.2)
X <- c(0.4, 0, 0.2, 0, 0.4)

test_that("tree.corr returns 0 without weights", {
  expect_equal(tree.corr.cond(P, Y, X, "Dn")[1], 0)
})

test_that("tree.corr returns expected values with consistent weights", {
  expect_equal(tree.corr.cond(P, Y, X, "Dn", weight = rep(1, 9))[1], 0.12, tolerance = 0.05)
  expect_equal(tree.corr.cond(P, Y, X, "Dn", weight = c(1, 2, 3, 4, 5, 6, 7, 8, 9))[1], 0.93, tolerance = 0.05)
})

test_that("tree.corr fails on invalid method", {
  expect_error(tree.corr.cond(P, Y, X, "invalid"))
})
