test_that("NodeDepth runs on a test tree", {
  num_tips <- 5
  
  tree <- ape::rtree(num_tips)
  alphaP <- c(rep(1, length(tree$tip.label)), rep(0, tree$Nnode))
  alphaQ <- c(rep(1, length(tree$tip.label)), rep(0, tree$Nnode))
  DataPQ <- DataGenerating(num_tips, num_tips, alphaP, alphaQ, 10)
  
  ND <- NodeDepth(tree)
  
  expect_equal(length(tree$tip.label) + tree$Nnode, length(ND))
})