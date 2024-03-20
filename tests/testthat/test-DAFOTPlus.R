test_that("IndDAFOT function works on highly correlated data", {
  d <- test_IndDAFOT
  step <- 49

  start_time <- Sys.time()
  res <- IndDAFOT(d$P, d$Y, d$tree)
  end_time <- Sys.time()

  print(paste((end_time - start_time) / step, " per step for IndDAFOT"))
  expect_true(all(res$P < 0.05))
})

test_that("ConIndDAFOT function works on highly correlated data", {
  d <- test_ConIndDAFOT
  step <- 49

  start_time_condgen <- Sys.time()
  res_condgen <- ConIndDAFOT(d$P, d$Y, d$X, d$tree, condgen = d$condgen, step = step)
  end_time_condgen <- Sys.time()

  print(paste((end_time_condgen - start_time_condgen) / step, " per step for ConIndDAFOT with condgen"))
  expect_true(all(res_condgen$P < 0.05))

  start_time_Ex <- Sys.time()
  res_Ex <- ConIndDAFOT(d$P, d$Y, d$X, d$tree, ExY = d$ExY, ExX = d$ExX, Exk = d$Exk, step = step)
  end_time_Ex <- Sys.time()

  print(paste((end_time_Ex - start_time_Ex) / step, " per step for ConIndDAFOT with Ex's"))
  expect_true(all(res_Ex$P < 0.05))
})
