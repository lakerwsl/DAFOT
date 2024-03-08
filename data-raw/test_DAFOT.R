set.seed(59049)

source("data-raw/util.R")

gg_db <- setup_gg_db()
subtree <- create_subtree(gg_db$taxa_dbi, gg_db$tree)
tree <- subtree$tree
tree$tip.label <- 1:247

### Model groups ###

mP <- mQ <- m <- 200
n <- 1000
delta <- 0.5

alphaP <- c(rep(1, length(tree$tip.label)), rep(0, tree$Nnode))
alphaQ <- c(rep(1, length(tree$tip.label)), rep(0, tree$Nnode))
alphaQ[169] <- alphaQ[169] + delta
alphaQ[170] <- alphaQ[170] - delta
DataPQ <- DataGenerating(mP, mQ, alphaP, alphaQ, n)
P <- DataPQ$P
Q <- DataPQ$Q

test_DAFOT <- list(P = P, Q = Q, tree = tree)
usethis::use_data(test_DAFOT, overwrite = TRUE)
