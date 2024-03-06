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
#' @param Tree A phylo class. Phylogenetic tree.
#' @param times A integer. The times of permutation test.
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
#' Tree <- rtree(100)
#' alphaP <- c(rep(1, length(Tree$tip.label)), rep(0, Tree$Nnode))
#' alphaQ <- c(rep(1, length(Tree$tip.label)), rep(0, Tree$Nnode))
#' DataPQ <- DataGenerating(100, 100, alphaP, alphaQ, 1000)
#' n <- 50
#' DAFOT(DataPQ$P, DataPQ$Q, Tree, n, 0.05)
#'
#' @author Shulei Wang
DAFOT <- function(P, Q, Tree, times = 200, alpha = 0.05) {
  AccuProbMt <- function(rP, Tree) {
    TTedge <- Tree$edge
    m <- length(Tree$tip.label) + Tree$Nnode
    Nodenum <- 1:m
    TTedge <- Tree$edge
    P <- rP

    ND <- NodeDepth(Tree)

    maxDepth <- max(ND)
    for (k in 2:maxDepth)
    {
      Tnode <- Nodenum[ND == (maxDepth - k + 2)]
      Index <- (TTedge[, 2] %in% Tnode)
      Tedge <- matrix(TTedge[Index, ], ncol = 2)
      for (i in 1:sum(Index))
      {
        P[Tedge[i, 1], ] <- P[Tedge[i, 1], ] + P[Tedge[i, 2], ]
      }
    }
    return(P)
  }

  NodeDepth <- function(Tree) {
    m <- length(Tree$tip.label) + Tree$Nnode
    Nodenum <- 1:m

    Troot <- tidytree::rootnode(as_tibble(Tree))
    Troot <- Troot$node
    NodeDepth <- rep(0, m)
    NodeDepth[Troot] <- 1
    Depth <- 1
    TTedge <- Tree$edge
    while (sum(NodeDepth <= 0) > 0) {
      Tnode <- Nodenum[NodeDepth == Depth]
      Index <- (TTedge[, 1] %in% Tnode)
      Tedge <- matrix(TTedge[Index, ], ncol = 2)
      Depth <- Depth + 1
      for (i in 1:sum(Index))
      {
        NodeDepth[Tedge[i, 2]] <- Depth
      }
    }
    return(NodeDepth)
  }

  standEdge <- function(EdgeP, EdgeQ, nP, nQ) {
    meanP <- apply(EdgeP, 1, mean)
    sdP <- apply(EdgeP, 1, stats::sd)
    meanQ <- apply(EdgeQ, 1, mean)
    sdQ <- apply(EdgeQ, 1, stats::sd)
    Diffmean <- meanP - meanQ
    Diffsd <- sqrt(sdP^2 / nP + sdQ^2 / nQ)

    Diffstat <- Diffmean
    Diffstat[Diffsd > 0] <- Diffmean[Diffsd > 0] / Diffsd[Diffsd > 0]
    Diffstat[Diffsd <= 0] <- 0
    Diffstat <- abs(Diffstat)

    return(Diffstat)
  }

  stopifnot(is.matrix(P))
  stopifnot(is.matrix(Q))
  stopifnot(is(Tree, "phylo"))
  times <- as.integer(times)
  stopifnot(is.numeric(alpha))
  stopifnot(alpha < 1 && alpha > 0)
  stopifnot(nrow(P) == nrow(Q))
  stopifnot(nrow(P) == length(Tree$tip.label) + Tree$Nnode)
  stopifnot(all(apply(P, 2, sum) == 1))
  stopifnot(all(apply(Q, 2, sum) == 1))
  stopifnot(all(P >= 0))
  stopifnot(all(Q >= 0))

  EdgeP <- AccuProbMt(P, Tree)
  EdgeQ <- AccuProbMt(Q, Tree)
  nP <- ncol(EdgeP)
  nQ <- ncol(EdgeQ)

  Diffstat <- standEdge(EdgeP, EdgeQ, nP, nQ)
  Stat <- max(Diffstat)

  EdgePQ <- cbind(EdgeP, EdgeQ)
  SimulatedRe <- rep(0, times)
  for (i in 1:times)
  {
    sampleQ <- 1:(nP + nQ)
    sampleP <- sample(sampleQ, nP)
    sampleQ <- sampleQ[!(sampleQ %in% sampleP)]
    sampleEdgeP <- EdgePQ[, sampleP]
    sampleEdgeQ <- EdgePQ[, sampleQ]

    sDiffstat <- standEdge(sampleEdgeP, sampleEdgeQ, nP, nQ)
    SimulatedRe[i] <- max(sDiffstat)
  }
  PValue <- (1 + sum(SimulatedRe > Stat)) / (1 + times)
  StatThre <- stats::quantile(SimulatedRe, 1 - alpha)

  Nodes <- 1:length(Diffstat)
  Nodes <- Nodes[Diffstat > StatThre]
  Edges <- Tree$edge
  ActiveEdge <- Edges[Edges[, 2] %in% Nodes, ]
  Re <- list(Stat = Stat, P = PValue, Thre = StatThre, Active = ActiveEdge)
  class(Re) <- "dafot"

  return(Re)
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
#' @param Tree A phylo class. Phylogenetic tree
#' @param t a numeric. The proportion of sample from 1st group.
#'
#' @return \code{SCalculation} returns a value for effctive number of tree.
#'
#' @export
#'
#' @examples
#' library(ape)
#' Tree <- rtree(100)
#' alphaP <- c(rep(1, length(Tree$tip.label)), rep(0, Tree$Nnode))
#' alphaQ <- c(rep(1, length(Tree$tip.label)), rep(0, Tree$Nnode))
#' DataPQ <- DataGenerating(100, 100, alphaP, alphaQ, 1000)
#' sdP <- apply(DataPQ$P, 1, sd)
#' sdQ <- apply(DataPQ$Q, 1, sd)
#' SCalculation(sdP, sdQ, Tree, 100 / (100 + 100))
#'
#' @author Shulei Wang
SCalculation <- function(sP, sQ, Tree, t) {
  AccuProb <- function(rP, Tree) {
    TTedge <- Tree$edge
    m <- length(Tree$tip.label) + Tree$Nnode
    Nodenum <- 1:m
    TTedge <- Tree$edge
    P <- rP

    ND <- NodeDepth(Tree)

    maxDepth <- max(ND)
    for (k in 2:maxDepth)
    {
      Tnode <- Nodenum[ND == (maxDepth - k + 2)]
      Index <- (TTedge[, 2] %in% Tnode)
      Tedge <- matrix(TTedge[Index, ], ncol = 2)
      for (i in 1:sum(Index))
      {
        P[Tedge[i, 1]] <- P[Tedge[i, 1]] + P[Tedge[i, 2]]
      }
    }
    return(P)
  }

  NodeDepth <- function(Tree) {
    m <- length(Tree$tip.label) + Tree$Nnode
    Nodenum <- 1:m

    Troot <- tidytree::rootnode(tibble::as_tibble(Tree))
    Troot <- Troot$node
    NodeDepth <- rep(0, m)
    NodeDepth[Troot] <- 1
    Depth <- 1
    TTedge <- Tree$edge
    while (sum(NodeDepth <= 0) > 0) {
      Tnode <- Nodenum[NodeDepth == Depth]
      Index <- (TTedge[, 1] %in% Tnode)
      Tedge <- matrix(TTedge[Index, ], ncol = 2)
      Depth <- Depth + 1
      for (i in 1:sum(Index))
      {
        NodeDepth[Tedge[i, 2]] <- Depth
      }
    }
    return(NodeDepth)
  }

  stopifnot(is.vector(sP))
  stopifnot(is.vector(sQ))
  stopifnot(is(Tree, "phylo"))
  stopifnot(is.numeric(t))
  stopifnot(t <= 1 && t > 0)
  stopifnot(length(sP) == length(sQ))
  stopifnot(length(sP) == length(Tree$tip.label) + Tree$Nnode)
  stopifnot(all(sP >= 0))
  stopifnot(all(sQ >= 0))

  EdgeP <- AccuProb(sP, Tree)
  EdgeQ <- AccuProb(sQ, Tree)
  ASigma <- EdgeP * (1 - t) + EdgeQ * t

  J <- 20
  s <- length(ASigma)
  Nodenum <- 1:s
  NDepth <- NodeDepth(Tree)
  STpi <- 0

  for (j in 1:J)
  {
    TThre <- 4 / (2^j)
    TNList <- Nodenum[(ASigma < TThre) & (ASigma >= TThre / 2)]
    TNListND <- NDepth[TNList]

    while (length(TNList) > 0) {
      TempMinIndex <- which(TNListND == max(TNListND))[1]
      TAncesor <- c(tidytree::ancestor(tibble::as_tibble(Tree), TNList[TempMinIndex])$node, TNList[TempMinIndex])
      TAncesor <- TAncesor[(ASigma[TAncesor] < TThre) & (ASigma[TAncesor] >= TThre / 2)]
      STpi <- STpi + 1

      TempTNLIndex <- TNList %in% TAncesor
      TNListND <- TNListND[!TempTNLIndex]
      TNList <- TNList[!TempTNLIndex]
    }
  }

  return(STpi)
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

  return(list(P = trP, Q = trQ))
}
