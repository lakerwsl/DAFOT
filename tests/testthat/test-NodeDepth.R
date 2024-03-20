test_that("NodeDepth runs on a test tree", {
  num_tips <- 5
  
  tree <- ape::rtree(num_tips)
  
  ND <- NodeDepth(tree)
  
  expect_equal(length(tree$tip.label) + tree$Nnode, length(ND))
})