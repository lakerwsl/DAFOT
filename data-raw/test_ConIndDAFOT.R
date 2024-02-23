set.seed(59049)

source("data-raw/util.R")

gg_db <- setup_gg_db()
subtree <- create_subtree(gg_db$taxa_dbi, gg_db$tree)
tree <- subtree$tree
Clade1 <- subtree$Clade1
Clade2 <- subtree$Clade2

### Model effect ###

f_s <- function(x) {
  sum(sin(2 * pi * x))
}
f_m <- function(x) {
  mean(cos(2 * pi * x))
}
f_y <- function(y) {
  cos(pi * y / 3)
}

condgen <- function(X) {
  f_s <- function(x) {
    sum(sin(2 * pi * x))
  }

  Y <- rnorm(1, mean = f_s(X))

  Y
}

n <- 100
m <- 5000
beta <- 100
model <- 1

X <- matrix(runif(n * 5), nrow = n)
means <- apply(X, 1, f_s)
Y <- rnorm(n, mean = means)
P <- matrix(0, nrow = length(tree$tip.label) + tree$Nnode, ncol = n)

if (model == 1) {
  means <- apply(X, 1, f_m) * 0 + f_y(Y) * 1
} else if (model == 2) {
  means <- apply(X, 1, f_m) * 0.2 + f_y(Y) * 0.8
} else if (model == 3) {
  means <- apply(X, 1, f_m) * 0.4 + f_y(Y) * 0.6
} else if (model == 4) {
  means <- apply(X, 1, f_m) * 0.6 + f_y(Y) * 0.4
} else if (model == 5) {
  means <- apply(X, 1, f_m) * 0.8 + f_y(Y) * 0.2
}

for (i in 1:n) {
  alphaQ <- c(rep(1, length(tree$tip.label)), rep(0, tree$Nnode))
  alphaQ[Clade1] <- alphaQ[Clade1] + means[i]
  alphaQ[Clade2] <- alphaQ[Clade2] - means[i]

  P[, i] <- gtools::rdirichlet(1, alphaQ)
  P[, i] <- rmultinom(1, m, P[, i])
  P[, i] <- P[, i] / sum(P[, i])
}

ExX <- matrix(runif(n * 5 * beta), nrow = n * beta)
means <- apply(ExX, 1, f_s)
ExY <- rnorm(n * beta, mean = means)
Exk <- 30

### Save data ###

test_ConIndDAFOT <- list(P = P, Y = Y, X = X, tree = tree, condgen = condgen, ExY = ExY, ExX = ExX, Exk = Exk)
usethis::use_data(test_ConIndDAFOT, overwrite = TRUE)
