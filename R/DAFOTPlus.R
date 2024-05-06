#' @title Phylogenetic Independence Test with Rank Correlation (IndDAFOT)
#'
#' @description A phylogenetic independence test directly measures the association between the outcome of interest and the total microbial abundance in each lineage (subtree) by a rank correlation.
#'
#' @details This method is designed to capture possible nonlinear associations between the outcome of interest and microbial abundance. It aggregates the rank correlations for different lineages in two ways: the weighted sum and maximum.
#'
#' @importFrom independence hoeffding.D.test
#' @importFrom independence tau.star.test
#' @importFrom independence hoeffding.refined.test
#' @importFrom tidytree rootnode
#' @importFrom ggtree ggtree
#' @importFrom ape rtree
#'
#' @param P A numeric matrix or data frame. Compositional data. Row represents nodes on the phylogenetic tree. Column represents samples.
#' @param Y A vector. The interested variable.
#' @param tree \code{phylo} class. Phylogenetic tree.
#' @param method Dn or Rn or Tn. Dn is the Hoeffding's D test; Rn is the Blum-Kiefer-Rosenblatt's R; Tn is the Bergsma-Dassios-Yanaginoto's tau test.
#' @param step Permutation times.
#'
#' @return a \code{list} with components:
#' \describe{
#' \item{Stat}{a vector with length 2. The first result is the correlation coefficient using the weighted sum approach; the second result is the correlation coefficient using maximum approach.}
#' \item{P}{a vector with length 2. The first result is the obtained p-value using the weighted sum approach; the second result is the obtained p-value using maximum approach.}
#' }
#'
#' @export
#'
#' @examples
#' library(ape)
#' Tree <- rtree(100)
#' P <- rbind(
#'   matrix(1, nrow = length(Tree$tip.label), ncol = 100),
#'   matrix(0, nrow = Tree$Nnode, ncol = 100)
#' )
#' Y <- c(rep(1, 50), rep(2, 50))
#' n <- 1
#' IndDAFOT(P, Y, Tree, step = n)
IndDAFOT <- function(P, Y, tree, method = "Dn", step = 200) {
  EdgeP <- AccuProb(P, tree)
  EdgeL <- EdgeLExtract(tree)
  EdgeP <- EdgeP[EdgeL > 0, ]
  EdgeL <- EdgeL[EdgeL > 0]

  n <- ncol(EdgeP)
  t <- tree.corr(EdgeP, Y, method, weight = EdgeL)
  sample_matrix <- replicate(step, sample(1:n, replace = F))
  Y1 <- matrix(Y[sample_matrix], nrow = n)
  Phi <- t(apply(Y1, 2, function(Y) tree.corr(EdgeP, Y, method, weight = EdgeL)))
  Pvalue <- c(mean(c(Phi[, 1], t[1]) >= t[1]), mean(c(Phi[, 2], t[2]) >= t[2]))

  list(Stat = t, P = Pvalue)
}


#' @title Phylogenetic Conditional Independence Test with Rank Correlation (ConIndDAFOT)
#'
#' @description A phylogenetic independence test directly measures the association between the outcome of interest and the total microbial abundance with adjustment of confounding effects.
#'
#' @details This method is designed to control confounding effect in the basis of IndDAFOT. Similar to the idea of the nearest neighbor mehtod, we consider the nearest neighbor conditional rank correlation coefficient for each lineage. And conditional rank correlations for lineage are aggregated as the weighted sum and maximum to capture signals.
#'
#' @importFrom stats dist
#' @importFrom independence tau.star.test
#' @importFrom independence hoeffding.refined.test
#' @importFrom tidytree rootnode
#' @importFrom ggtree ggtree
#' @importFrom RANN nn2
#' @importFrom ape rtree
#'
#' @param P A numeric matrix or data frame. Compositional data. Row represents nodes on the phylogenetic tree. Column represents samples.
#' @param Y A vector. The interested variable.
#' @param X The confounding variables. Row represents samples. Each column corresponds to each confounding variable.
#' @param Tree \code{phylo} class. Phylogenetic tree.
#' @param condgen \code{function}. Conditional distribution generating function.
#' @param ExY Extra data generated from the conditional distribution of Y given X.
#' @param ExX Extra data generated from the distribution of X.
#' @param Exk the number of neighbors to draw Y from the neighbors of X.
#' @param method Dn or Rn or Tn. Dn is the Hoeffding's D test; Rn is the Blum-Kiefer-Rosenblatt's R; Tn is the Bergsma-Dassios-Yanaginoto's tau test.
#' @param neighbor The number of neighbors to estimate local rank correlation.
#' @param step Permutation times.

#'
#' @return a \code{list} with components:
#' \describe{
#' \item{Stat}{a vector with length 2. The first result is the correlation coefficient using the weighted sum approach; the second result is the correlation coefficient using maximum approach.}
#' \item{P}{a vector with length 2. The first result is the obtained p-value using the weighted sum approache; the second result is the obtained p-value using maximum approach.}
#' }
#'
#' @export
#'
#' @examples
#' library(ape)
#' Tree <- rtree(100)
#' P <- rbind(
#'   matrix(1, nrow = length(Tree$tip.label), ncol = 100),
#'   matrix(0, nrow = Tree$Nnode, ncol = 100)
#' )
#' Y <- c(rep(1, 50), rep(2, 50))
#' X <- matrix(rnorm(500), nrow = 100)
#' n <- 1
#' ConIndDAFOT(P, Y, X, Tree, step = n)
ConIndDAFOT <- function(P, Y, X, Tree, condgen = NULL, ExY = NULL, ExX = NULL, Exk = min(10, length(ExY)), method = "Dn", neighbor = min(10, length(Y)), step = 200) {
  EdgeP <- AccuProb(P, Tree)
  EdgeL <- EdgeLExtract(Tree)
  EdgeP <- EdgeP[EdgeL > 0, ]
  EdgeL <- EdgeL[EdgeL > 0]

  n <- ncol(EdgeP)
  t <- tree.corr.cond(EdgeP, Y, X, method, weight = EdgeL, neighbor = neighbor)

  Phi <- matrix(0, nrow = step, ncol = 2)

  if (!is.null(condgen)) {
    for (i in 1:step) {
      NewY <- sapply(1:nrow(X), function(i, X) {
        condgen(X[i, ])
      }, X)
      Phi[i, ] <- tree.corr.cond(EdgeP, NewY, X, method, weight = EdgeL, neighbor = neighbor)
    }
  } else if (!is.null(ExY) & !is.null(ExX)) {
    ExNN <- RANN::nn2(ExX, X, k = Exk)$nn.idx

    for (i in 1:step) {
      NewIndex <- apply(ExNN, 1, sample, 1)
      NewY <- ExY[NewIndex]
      Phi[i, ] <- tree.corr.cond(EdgeP, NewY, X, method, weight = EdgeL, neighbor = neighbor)
    }
  } else {
    distm <- as.matrix(stats::dist(X))
    NewIndex <- sapply(1:ncol(distm), function(X) order(distm[, X], decreasing = F)[1:neighbor])

    for (i in 1:step) {
      Ind <- apply(NewIndex, 2, function(x) sample(x, 1))
      NewY <- Y[Ind]
      Phi[i, ] <- tree.corr.cond(EdgeP, NewY, X, method, weight = EdgeL, neighbor = neighbor)
    }
  }

  PValue <- c(mean(c(Phi[, 1], t[1]) >= t[1]), mean(c(Phi[, 2], t[2]) >= t[2]))

  list(Stat = t, P = PValue)
}
