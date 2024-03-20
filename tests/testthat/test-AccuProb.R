test_that("AccuProb converts node abundance to edge abundance correctly", {
  # Create a simple phylogenetic tree
  tree <- ggtree::read.tree(text = "(A:1,B:1,(C:1,D:1):1);")

  # Create a sample data frame representing node abundance
  node_abundance <- data.frame(A = c(1, 2), B = c(3, 4), C = c(5, 6), D = c(7, 8))

  # Call the AccuProb function to convert node abundance to edge abundance
  edge_abundance <- AccuProb(node_abundance, tree)

  # Expected edge abundance
  expected_edge_abundance <- data.frame(A = c(3, 3), C = c(5, 5), D = c(7, 7))

  # Assert that the result matches the expected edge abundance
  expect_equal(edge_abundance, expected_edge_abundance)
})
