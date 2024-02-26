set.seed(59049)

source("data-raw/util.R")

gg_db <- setup_gg_db()
subtree <- create_subtree(gg_db$taxa_dbi, gg_db$tree)
tree <- subtree$tree
Clade1 <- subtree$Clade1
Clade2 <- subtree$Clade2

### Model groups ###

n <- 200
m <- 5000

P <- matrix(0, nrow = length(tree$tip.label) + tree$Nnode, ncol = n)
Q <- matrix(0, nrow = length(tree$tip.label) + tree$Nnode, ncol = n)

for (i in 1:n) {
  alphaP <- alphaQ <- c(rep(1, length(tree$tip.label)), rep(1, tree$Nnode))
  alphaP[Clade1] <- alphaP[Clade1] + runif(length(alphaP[Clade1]), max = 5)
  alphaQ[Clade1] <- alphaQ[Clade1] + runif(length(alphaQ[Clade1]), max = 5)

  P[, i] <- rdirichlet(1, alphaP)
  Q[, i] <- rdirichlet(1, alphaQ)
  P[, i] <- rmultinom(1, m, P[, i])
  Q[, i] <- rmultinom(1, m, Q[, i])
  P[, i] <- P[, i] / sum(P[, i])
  Q[, i] <- Q[, i] / sum(Q[, i])
}

test_DAFOT <- list(P = P, Q = Q, tree = tree)
usethis::use_data(test_DAFOT, overwrite = TRUE)
