test_that("IndDAFOT function works as expected", {
  d <- test_IndDAFOT
  res <- IndDAFOT(d$P, d$Y, d$tree)

  expect_true(all(res$P < 0.05))
})

# ~5s per iteration (step value) on test data
test_that("ConIndDAFOT function works as expected", {
  d <- test_ConIndDAFOT
  # res_condgen <- ConIndDAFOT(d$P, d$Y, d$X, d$tree, condgen = d$condgen, step = 5)



  res_Ex <- ConIndDAFOT(d$P, d$Y, d$X, d$tree, ExY = d$ExY, ExX = d$ExX, Exk = d$Exk, step = 5)

  expect_true(all(res_Ex$P < 0.05))
})
