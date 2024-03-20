test_that("AccuProb runs on a random test tree", {
  num_tips <- 25

  tree <- ape::rtree(num_tips)

  node_abundance <- data.frame(matrix(rnorm(num_tips), nrow = 1))

  edge_abundance <- AccuProb(node_abundance, tree)

  expect_equal(nrow(edge_abundance), tree$Nnode + length(tree$tip.label))
})

test_that("AccuProb returns the expected answer for a known tree and node abundance data frame", {
  tree <- structure(list(edge = structure(c(6L, 7L, 7L, 6L, 8L, 9L, 9L, 8L, 7L, 1L, 2L, 8L, 9L, 3L, 4L, 5L), dim = c(8L, 2L)), tip.label = c("t5", "t4", "t2", "t1", "t3"), Nnode = 4L, edge.length = c(0.759698169771582, 0.237785059725866, 0.974591086851433, 0.0666242216248065, 0.423836079193279, 0.51583820884116, 0.532561222324148, 0.403174153761938)), class = "phylo", order = "cladewise")

  node_abundance <- structure(c(0.4, 0, 0.4, 0, 0.2, 0, 0, 0, 0, 0, 0.4, 0.2, 0.4, 0, 0, 0, 0, 0, 0.8, 0, 0, 0.2, 0, 0, 0, 0, 0, 0.8, 0, 0, 0, 0.2, 0, 0, 0, 0, 0.4, 0.4, 0, 0.2, 0, 0, 0, 0, 0), dim = c(9L, 5L))

  edge_abundance <- AccuProb(node_abundance, tree)

  expect_equal(edge_abundance, structure(c(0.4, 0, 0.4, 0, 0.2, 1, 0.4, 0.6, 0.4, 0, 0.4, 0.2, 0.4, 0, 1, 0.4, 0.6, 0.6, 0.8, 0, 0, 0.2, 0, 1, 0.8, 0.2, 0.2, 0.8, 0, 0, 0, 0.2, 1, 0.8, 0.2, 0, 0.4, 0.4, 0, 0.2, 0, 1, 0.8, 0.2, 0.2), dim = c(9L, 5L)))
})
