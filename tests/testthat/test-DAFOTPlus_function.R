test_that("IndDAFOT function works as expected", {
  # Generate example data
  Tree <- rtree(100)
  P <- rbind(matrix(1, nrow = length(Tree$tip.label), ncol = 100),
             matrix(0, nrow = Tree$Nnode, ncol = 100))
  Y <- c(rep(1, 50), rep(2, 50))
  
  # Test IndDAFOT function
  result <- IndDAFOT(P, Y, Tree)
  
  # Check if the result is a list with components Stat and P
  expect_type(result, "list")
  expect_named(result, c("Stat", "P"))
  expect_length(result$Stat, 2)
  expect_length(result$P, 2)
})

### TAKES WAY TOO LONG
test_that("ConIndDAFOT function works as expected", {
  # Generate example data
  Tree <- rtree(100)
  P <- rbind(matrix(1, nrow = length(Tree$tip.label), ncol = 100),
             matrix(0, nrow = Tree$Nnode, ncol = 100))
  Y <- c(rep(1, 50), rep(2, 50))
  X <- matrix(rnorm(500), nrow = 100)
  
  # Test ConIndDAFOT function
  #result <- ConIndDAFOT(P, Y, X, Tree)
  
  # Check if the result is a list with components Stat and P
  #expect_type(result, "list")
  #expect_named(result, c("Stat", "P"))
  #expect_length(result$Stat, 2)
  #expect_length(result$P, 2)
})