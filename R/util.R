#' Return the depth of each node on the tree
#'
#' @importFrom ape nodepath
#'
#' @param Tree \code{phylo} class. Phylogenetic tree.
#'
#' @return An integer vector containing the depth of each node.
NodeDepth <- function(Tree) {
  np <- nodepath(Tree)
  depths <- rep(0, length(Tree$tip.label) + Tree$Nnode)
  for (i in 1:length(np)) {
    for (j in 1:length(np[[i]])) {
      depths[np[[i]][j]] <- j
    }
  }

  depths
}

#' Convert node abundance to edge abundance
#'
#' @param P A numeric matrix or data frame. Compositional data. Row represents nodes on the phylogenetic tree. Columns represents samples.
#' @param Tree \code{phylo} class. Phylogenetic tree.
#' @param r If calculation is being performed on variance or values.
#'
#' @return A numeric matrix or data frame. Compositional data. Row represents edges on the phylogenetic tree. Columns represent samples.
AccuProb <- function(P, Tree, r = FALSE) {
  TTedge <- Tree$edge
  m <- length(Tree$tip.label) + Tree$Nnode
  Nodenum <- 1:m
  ND <- NodeDepth(Tree)
  maxDepth <- max(ND)

  for (k in 2:maxDepth)
  {
    Tnode <- Nodenum[ND == (maxDepth - k + 2)]
    Index <- (TTedge[, 2] %in% Tnode)
    if (r) {
      Tedge <- matrix(TTedge[Index, ], ncol = 2)

      for (i in 1:sum(Index))
      {
        P[Tedge[i, 1]] <- P[Tedge[i, 1]] + P[Tedge[i, 2]]
      }
    } else {
      Tedge <- TTedge[Index, , drop = FALSE]

      for (i in 1:sum(Index))
      {
        P[Tedge[i, 1], ] <- P[Tedge[i, 1], ] + P[Tedge[i, 2], ]
      }
    }
  }

  P
}

#' Extract sorted edge weights from the phylogenetic tree
#' Add zero for the root node edge length
#' Return all 0's if tree is invalid
#'
#' @importFrom tidytree rootnode
#'
#' @param Tree \code{phylo} class. Phylogenetic tree.
#'
#' @return A vector of numerics representing edge weights on the tree.
EdgeLExtract <- function(Tree) {
  Trank <- sort(Tree$edge[, 2], index.return = TRUE)
  EdgeL <- Tree$edge.length[Trank$ix]
  rootnum <- rootnode(Tree)

  if (rootnum > length(EdgeL)) {
    EdgeL <- c(EdgeL, 0)
  } else {
    EdgeL <- c(EdgeL[1:(rootnum - 1)], 0, EdgeL[rootnum:length(EdgeL)])
  }

  EdgeL
}


#' Standardize edge weights
#'
#' @importFrom stats sd
#'
#' @param EdgeP A numeric matrix or data frame. Compositional data. Row represents edges on the phylogenetic tree. Columns represent samples.
#' @param EdgeQ A numeric matrix or data frame. Compositional data. Row represents edges on the phylogenetic tree. Columns represent samples.
#' @param nP The number of columns in EdgeP
#' @param nQ The number of columns in EdgeQ
#'
#' @return A numeric vector. The differences between input edge matrixes.
standEdge <- function(EdgeP, EdgeQ, nP, nQ) {
  meanP <- apply(EdgeP, 1, mean)
  sdP <- apply(EdgeP, 1, sd)
  meanQ <- apply(EdgeQ, 1, mean)
  sdQ <- apply(EdgeQ, 1, sd)
  Diffmean <- meanP - meanQ
  Diffsd <- sqrt(sdP^2 / nP + sdQ^2 / nQ)

  Diffstat <- Diffmean
  Diffstat[Diffsd > 0] <- Diffmean[Diffsd > 0] / Diffsd[Diffsd > 0]
  Diffstat[Diffsd <= 0] <- 0
  Diffstat <- abs(Diffstat)

  Diffstat
}

#' Correlation between one single edge and Y
#'
#' @importFrom independence hoeffding.D.test tau.star.test hoeffding.refined.test
#'
#' @param EdgeP A numeric matrix or data frame. Compositional data. Row represents edges on the phylogenetic tree. Columns represent samples.
#' @param Y A vector. The interested variable.
#' @param method Dn or Rn or Tn. Dn is the Hoeffding's D test; Rn is the Blum-Kiefer-Rosenblatt's R; Tn is the Bergsma-Dassios-Yanaginoto's tau test.
#' @param weight A vector. Weights to apply to the test results.
#'
#' @return A vector. First field is the weighted sum of test results. Second is the max result.
tree.corr <- function(EdgeP, Y, method, weight = F) {
  n <- ncol(EdgeP)

  if (method == "Dn") {
    phi <- apply(EdgeP, 1, function(P) hoeffding.D.test(P, Y, precision = 1, collisions = F)$Dn)
  } else if (method == "Tn") {
    phi <- apply(EdgeP, 1, function(P) tau.star.test(P, Y, precision = 1, collisions = F)$Tn)
  } else if (method == "Rn") {
    phi <- apply(EdgeP, 1, function(P) hoeffding.refined.test(P, Y, precision = 1, collisions = F)$Rn)
  } else {
    stop("Please choose valid method.")
  }

  c(sum(weight * phi), max(phi))
}

#' Conditional correlation between one single edge and Y
#'
#' @importFrom independence hoeffding.D.test tau.star.test hoeffding.refined.test
#' @importFrom stats dist
#'
#' @param EdgeP A numeric matrix or data frame. Compositional data. Row represents edges on the phylogenetic tree. Columns represent samples.
#' @param Y A vector. The interested variable.
#' @param X The confounding variables. Row represents samples. Each column corresponds to each confounding variable.
#' @param method Dn or Rn or Tn. Dn is the Hoeffding's D test; Rn is the Blum-Kiefer-Rosenblatt's R; Tn is the Bergsma-Dassios-Yanaginoto's tau test.
#' @param neighbor The number of neighbors to estimate local rank correlation.
#' @param weight A vector. Weights to apply to the test results.
#'
#' @return A vector. First field is the weighted sum of test results. Second is the max result.
tree.corr.cond <- function(EdgeP, Y, X, method, neighbor = min(10, length(Y)), weight = F) {
  n <- ncol(EdgeP)
  d <- nrow(EdgeP)
  distm <- as.matrix(dist(X))
  phi <- rep(0, d)
  i.list <- apply(distm, 2, function(X) order(X, decreasing = F)[1:neighbor])
  Y1 <- matrix(Y[i.list], nrow = nrow(i.list))

  if (method == "Dn") {
    phi <- rowSums(apply(i.list, 2, function(x) apply(EdgeP[, x], 1, function(P) hoeffding.D.test(P, Y[x], precision = 1, collisions = F)$Dn)))
  } else if (method == "Tn") {
    phi <- rowSums(apply(i.list, 2, function(x) apply(EdgeP[, x], 1, function(P) tau.star.test(P, Y[x], precision = 1, collisions = F)$Tn)))
  } else if (method == "Rn") {
    phi <- rowSums(apply(i.list, 2, function(x) apply(EdgeP[, x], 1, function(P) hoeffding.refined.test(P, Y[x], precision = 1, collisions = F)$Rn)))
  } else {
    stop("Please choose valid method.")
  }

  phi <- phi / n

  c(sum(weight * phi), max(phi))
}
