test_that("NodeDepth runs on a random test tree", {
  num_tips <- 25

  tree <- ape::rtree(num_tips)

  ND <- NodeDepth(tree)

  expect_equal(length(tree$tip.label) + tree$Nnode, length(ND))
})

test_that("NodeDepth returns the expected answer for a known tree", {
  tree <- structure(list(edge = structure(c(6L, 7L, 7L, 6L, 8L, 9L, 9L, 8L, 7L, 1L, 2L, 8L, 9L, 3L, 4L, 5L), dim = c(8L, 2L)), tip.label = c("t5", "t4", "t2", "t1", "t3"), Nnode = 4L, edge.length = c(0.759698169771582, 0.237785059725866, 0.974591086851433, 0.0666242216248065, 0.423836079193279, 0.51583820884116, 0.532561222324148, 0.403174153761938)), class = "phylo", order = "cladewise")

  ND <- NodeDepth(tree)

  expect_equal(ND, c(3, 3, 4, 4, 3, 1, 2, 2, 3))
})
