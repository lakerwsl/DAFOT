# Unit tests for DAFOT function
test_that("DAFOT function returns expected results", {
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
  expect_is(dafot_result, "dafot")
  
  # Check if all components are present
  expect_named(dafot_result, c("Stat", "P", "Thre", "Active"))
  
  # Check if the Stat component is numeric
  expect_is(dafot_result$Stat, "numeric")
  
  # Check if the P component is numeric
  expect_is(dafot_result$P, "numeric")
  
  # Check if the Thre component is numeric
  expect_is(dafot_result$Thre, "numeric")
  
  # Check if the Active component is a matrix
  expect_is(dafot_result$Active, "matrix")
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
  s_calculation_result <- SCalculation(sdP, sdQ, Tree, 100/(100 + 100))
  
  # Check if the result is numeric
  expect_is(s_calculation_result, "numeric")
})

# Unit tests for DataGenerating function
test_that("DataGenerating function returns expected results", {
  # Test DataGenerating function
  data_generated <- DataGenerating(100, 100, rep(100, 1), rep(100, 2), 10000)
  
  # Check if the result is a list
  expect_is(data_generated, "list")
  
  # Check if the list contains two matrices
  expect_length(data_generated, 2)
  expect_is(data_generated$P, "matrix")
  expect_is(data_generated$Q, "matrix")
})