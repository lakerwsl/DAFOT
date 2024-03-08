set.seed(59049)

source("data-raw/util.R")

gg_db <- setup_gg_db()
subtree <- create_subtree(gg_db$taxa_dbi, gg_db$tree)
tree <- subtree$tree
Clade1 <- subtree$Clade1
Clade2 <- subtree$Clade2

### Model effect ###

delta <- 0.8
m <- 5000
n <- 200

Y <- runif(n)
P <- matrix(0, nrow = length(tree$tip.label) + tree$Nnode, ncol = n)
for (i in 1:n) {
  alphaQ <- c(rep(1, length(tree$tip.label)), rep(0, tree$Nnode))
  if (model == 1) {
    alphaQ[Clade1] <- alphaQ[Clade1] + delta * Y[i]
    alphaQ[Clade2] <- alphaQ[Clade2] - delta * Y[i]
  } else if (model == 2) {
    alphaQ[Clade1] <- alphaQ[Clade1] + delta * sin(Y[i] * 2 * pi)
    alphaQ[Clade2] <- alphaQ[Clade2] - delta * sin(Y[i] * 2 * pi)
  } else if (model == 3) {
    alphaQ[Clade1] <- alphaQ[Clade1] + delta * sin(Y[i] * 4 * pi)
    alphaQ[Clade2] <- alphaQ[Clade2] - delta * sin(Y[i] * 4 * pi)
  }
  P[, i] <- rdirichlet(1, alphaQ)
  P[, i] <- rmultinom(1, m, P[, i])
  P[, i] <- P[, i] / sum(P[, i])
}

### Save data ###

test_IndDAFOT <- list(P = P, Y = Y, tree = tree)
usethis::use_data(test_IndDAFOT, overwrite = TRUE)
