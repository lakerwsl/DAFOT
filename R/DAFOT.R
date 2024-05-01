#' @title Detector of Active Flow on a Tree
#'
#' @description \code{DAFOT} is a funtion used to detect the difference between two group of compositional data.
#'
#' @details
#' This function is developed based on an optimal transport perspective. It can help detect the specific location of difference at a tree.
#'
#' @seealso \code{\link{SCalculation}}
#'
#' @importFrom tidytree rootnode
#' @importFrom stats quantile
#' @importFrom stats sd
#' @importFrom methods is
#' @importFrom tibble as_tibble
#' @importFrom ape rtree
#'
#' @param P A matrix. Compositional data of 1st group, each row corresponds to each node of tree and each column corresponds to each sample.
#' @param Q A matrix. Compositional data of 2nd group, each row corresponds to each node of tree and each column corresponds to each sample.
#' @param tree A phylo class. Phylogenetic tree.
#' @param step A integer. The step of permutation test.
#' @param alpha A numeric. Significance level
#'
#' @return \code{DAFOT} returns an object of class \code{dafot}.
#'
#' An object of class \code{dafot} is a list containing following components:
#' \itemize{
#' \item \code{Stat} the maximum of test statistics
#' \item \code{P} P value calculated from permutation test
#' \item \code{Thre} the threshold for alpha level test
#' \item \code{Active} the set of active edges
#' }
#'
#' @export
#'
#' @examples
#' library(ape)
#' tree <- rtree(100)
#' alphaP <- c(rep(1, length(tree$tip.label)), rep(0, tree$Nnode))
#' alphaQ <- c(rep(1, length(tree$tip.label)), rep(0, tree$Nnode))
#' DataPQ <- DataGenerating(100, 100, alphaP, alphaQ, 1000)
#' n <- 50
#' DAFOT(DataPQ$P, DataPQ$Q, tree, n, 0.05)
#'
#' @author Shulei Wang
DAFOT <- function(P, Q, tree, step = 200, alpha = 0.05) {
  stopifnot(is.matrix(P))
  stopifnot(is.matrix(Q))
  stopifnot(is(tree, "phylo"))
  step <- as.integer(step)
  stopifnot(is.numeric(alpha))
  stopifnot(alpha < 1 && alpha > 0)
  stopifnot(nrow(P) == nrow(Q))
  stopifnot(nrow(P) == length(tree$tip.label) + tree$Nnode)
  print(apply(P, 2, sum))
  dput(apply(P, 2, sum))
  print(var(apply(P, 2, sum)))
  print(all(apply(P, 2, sum) == 1))
  stopifnot(var(apply(P, 2, sum)) == 0)
  stopifnot(var(apply(Q, 2, sum)) == 0)
  #stopifnot(all(apply(P, 2, sum) == 1))
  #stopifnot(all(apply(Q, 2, sum) == 1))
  stopifnot(all(P >= 0))
  stopifnot(all(Q >= 0))

  EdgeP <- AccuProb(P, tree)
  EdgeQ <- AccuProb(Q, tree)
  nP <- ncol(EdgeP)
  nQ <- ncol(EdgeQ)

  Diffstat <- standEdge(EdgeP, EdgeQ, nP, nQ)
  Stat <- max(Diffstat)

  EdgePQ <- cbind(EdgeP, EdgeQ)
  SimulatedRe <- rep(0, step)

  for (i in 1:step)
  {
    sampleQ <- 1:(nP + nQ)
    sampleP <- sample(sampleQ, nP)
    sampleQ <- sampleQ[!(sampleQ %in% sampleP)]
    sampleEdgeP <- EdgePQ[, sampleP]
    sampleEdgeQ <- EdgePQ[, sampleQ]

    sDiffstat <- standEdge(sampleEdgeP, sampleEdgeQ, nP, nQ)
    SimulatedRe[i] <- max(sDiffstat)
  }

  PValue <- (1 + sum(SimulatedRe > Stat)) / (1 + step)
  StatThre <- stats::quantile(SimulatedRe, 1 - alpha)

  Nodes <- 1:length(Diffstat)
  Nodes <- Nodes[Diffstat > StatThre]
  Edges <- tree$edge
  ActiveEdge <- Edges[Edges[, 2] %in% Nodes, ]
  Re <- list(Stat = Stat, P = PValue, Thre = StatThre, Active = ActiveEdge)
  class(Re) <- "dafot"

  Re
}


#' @title The effctive number of tree calculation
#'
#' @description \code{SCalculation} is used to calculate the effctive number of tree.
#'
#' @details
#' This function is used to calculate the effctive number of tree.
#'
#' @seealso \code{\link{DAFOT}}
#'
#' @importFrom tidytree rootnode
#' @importFrom tidytree ancestor
#' @importFrom methods is
#' @importFrom tibble as_tibble
#'
#' @param sP a vector. The component-wise variance of 1st group, each entry correpsonds each node of tree.
#' @param sQ a vector. The component-wise variance of 2nd group, each entry correpsonds each node of tree.
#' @param tree A phylo class. Phylogenetic tree
#' @param t a numeric. The proportion of sample from 1st group.
#'
#' @return \code{SCalculation} returns a value for effctive number of tree.
#'
#' @export
#'
#' @examples
#' library(ape)
#' tree <- rtree(100)
#' alphaP <- c(rep(1, length(tree$tip.label)), rep(0, tree$Nnode))
#' alphaQ <- c(rep(1, length(tree$tip.label)), rep(0, tree$Nnode))
#' DataPQ <- DataGenerating(100, 100, alphaP, alphaQ, 1000)
#' sdP <- apply(DataPQ$P, 1, sd)
#' sdQ <- apply(DataPQ$Q, 1, sd)
#' SCalculation(sdP, sdQ, tree, 100 / (100 + 100))
#'
#' @author Shulei Wang
SCalculation <- function(sP, sQ, tree, t) {
  stopifnot(is.vector(sP))
  stopifnot(is.vector(sQ))
  stopifnot(is(tree, "phylo"))
  stopifnot(is.numeric(t))
  stopifnot(t <= 1 && t > 0)
  stopifnot(length(sP) == length(sQ))
  stopifnot(length(sP) == length(tree$tip.label) + tree$Nnode)
  stopifnot(all(sP >= 0))
  stopifnot(all(sQ >= 0))

  EdgeP <- AccuProb(sP, tree, TRUE)
  EdgeQ <- AccuProb(sQ, tree, TRUE)
  ASigma <- EdgeP * (1 - t) + EdgeQ * t

  J <- 20
  s <- length(ASigma)
  Nodenum <- 1:s
  NDepth <- NodeDepth(tree)
  STpi <- 0

  for (j in 1:J)
  {
    TThre <- 4 / (2^j)
    TNList <- Nodenum[(ASigma < TThre) & (ASigma >= TThre / 2)]
    TNListND <- NDepth[TNList]

    while (length(TNList) > 0) {
      TempMinIndex <- which(TNListND == max(TNListND))[1]
      TAncesor <- c(tidytree::ancestor(tibble::as_tibble(tree), TNList[TempMinIndex])$node, TNList[TempMinIndex])
      TAncesor <- TAncesor[(ASigma[TAncesor] < TThre) & (ASigma[TAncesor] >= TThre / 2)]
      STpi <- STpi + 1

      TempTNLIndex <- TNList %in% TAncesor
      TNListND <- TNListND[!TempTNLIndex]
      TNList <- TNList[!TempTNLIndex]
    }
  }

  STpi
}


#' @title Generate random data on tree
#'
#' @description \code{DataGenerating} generate data from Dirichlet-multinomial distribution.
#'
#' @details
#' This function is used to generate m distributions from Dirichlet, and draw samples from multinomial distribution.
#'
#' @seealso \code{\link{DAFOT}}
#'
#' @importFrom gtools rdirichlet
#' @importFrom stats rmultinom
#'
#' @param mP an integer. The number of samples in the first group
#' @param mQ an integer. The number of samples in the second group
#' @param alphaP A vector of numeric. the parameter in Dirichlet distribution for the first group, each entry corresponds to each node of tree.
#' @param alphaQ A vector of numeric. the parameter in Dirichlet distribution for the second group, each entry corresponds to each node of tree.
#' @param n an integer. The number of reads drawn for each sample.
#'
#' @return \code{DataGenerating} returns a list of two matrix of the raw data generated from Dirichlet-multinomial distribution.
#'
#' @export
#'
#' @examples
#' alphaP <- rep(100, 1)
#' alphaQ <- rep(100, 2)
#' DataGenerating(100, 100, alphaP, alphaQ, 10000)
#'
#' @author Shulei Wang
DataGenerating <- function(mP, mQ, alphaP, alphaQ, n) {
  trP <- t(gtools::rdirichlet(mP, alphaP))
  trQ <- t(gtools::rdirichlet(mQ, alphaQ))

  for (i in 1:mP)
  {
    trP[, i] <- stats::rmultinom(1, n, trP[, i])
    trP[, i] <- trP[, i] / sum(trP[, i])
  }

  for (i in 1:mQ)
  {
    trQ[, i] <- stats::rmultinom(1, n, trQ[, i])
    trQ[, i] <- trQ[, i] / sum(trQ[, i])
  }

  list(P = trP, Q = trQ)
}
