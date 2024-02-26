test_that("DAFOT function works on highly correlated data", {
  d <- test_DAFOT
  step <- 199

  start_time <- Sys.time()
  res <- DAFOT(d$P, d$Q, d$tree, step, 0.05)
  end_time <- Sys.time()

  print(paste((end_time - start_time) / step, " per step for DAFOT"))
  # expect_true(res$P == 1 / (1 + step))
})

test_that("DAFOT function returns poorly correlated results on randomly generated data", {
  # Generate example data
  Tree <- rtree(100)
  alphaP <- c(rep(1, length(Tree$tip.label)), rep(0, Tree$Nnode))
  alphaQ <- c(rep(1, length(Tree$tip.label)), rep(0, Tree$Nnode))
  DataPQ <- DataGenerating(100, 100, alphaP, alphaQ, 1000)
  P <- DataPQ$P
  Q <- DataPQ$Q

  # Test DAFOT function
  dafot_result <- DAFOT(P, Q, Tree, 100, 0.05)

  # Check the class of the result
  expect_s3_class(dafot_result, "dafot")

  # Check if all components are present
  expect_named(dafot_result, c("Stat", "P", "Thre", "Active"))

  # Check Stat component
  expect_type(dafot_result$Stat, "double")
  expect_gt(dafot_result$Stat, 1)

  # Check P component
  expect_type(dafot_result$P, "double")
  expect_gt(dafot_result$P, 0.05)

  # Check Thre component
  expect_type(dafot_result$Thre, "double")
})

# Unit tests for SCalculation function
test_that("SCalculation function returns expected results", {
  # Generate example data
  Tree <- rtree(100)
  alphaP <- c(rep(1, length(Tree$tip.label)), rep(0, Tree$Nnode))
  alphaQ <- c(rep(1, length(Tree$tip.label)), rep(0, Tree$Nnode))
  DataPQ <- DataGenerating(100, 100, alphaP, alphaQ, 1000)
  sdP <- apply(DataPQ$P, 1, sd)
  sdQ <- apply(DataPQ$Q, 1, sd)

  # Test SCalculation function
  s_calculation_result <- SCalculation(sdP, sdQ, Tree, 100 / (100 + 100))

  # Check if the result is numeric
  expect_type(s_calculation_result, "double")
})

# Unit tests for DataGenerating function
test_that("DataGenerating function generates expected data structure", {
  # Generate example data
  mP <- 100
  mQ <- 100
  alphaP <- rep(100, 1)
  alphaQ <- rep(100, 2)
  n <- 10000

  # Test DataGenerating function
  data_generated <- DataGenerating(mP, mQ, alphaP, alphaQ, n)

  # Check if the result is a list
  expect_type(data_generated, "list")

  # Check if the dimensions of matrices are correct
  expect_equal(dim(data_generated$P), c(length(alphaP), mP))
  expect_equal(dim(data_generated$Q), c(length(alphaQ), mQ))

  # Check if probabilities sum to 1 for each sample
  expect_true(all(abs(apply(data_generated$P, 2, sum) - 1) < 1e-6))
  expect_true(all(abs(apply(data_generated$Q, 2, sum) - 1) < 1e-6))
})
