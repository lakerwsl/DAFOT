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

# Return the depth of each node on the tree.
NodeDepth <- function(Tree) {
  m <- length(Tree$tip.label) + Tree$Nnode
  Nodenum <- 1:m
  Troot <- tidytree::rootnode(Tree)
  NodeDepth <- rep(0, m)
  NodeDepth[Troot] <- 1
  Depth <- 1
  TTedge <- Tree$edge
  
  while (sum(NodeDepth <= 0) > 0) {
    Tnode <- Nodenum[NodeDepth == Depth]
    Index <- (TTedge[, 1] %in% Tnode)
    
    if (sum(Index) > 0) {
      Tedge <- TTedge[Index, , drop = FALSE]
      Depth <- Depth + 1
      NodeDepth[Tedge[, 2]] <- Depth
    }
  }
  
  NodeDepth
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
  
  NodeDepth
}

# Extract edge weights from the phylogenetic tree.
EdgeLExtract <- function(Tree) {
  Trank <- sort(Tree$edge[, 2], index.return = TRUE)
  EdgeL <- Tree$edge.length[Trank$ix]
  rootnum <- tidytree::rootnode(Tree)
  
  if (rootnum > length(EdgeL)) {
    EdgeL <- c(EdgeL, 0)
  } else {
    EdgeL <- c(EdgeL[1:(rootnum - 1)], 0, EdgeL[rootnum:length(EdgeL)])
  }
  
  EdgeL
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
  
  Diffstat
}

# Correlation between one single edge and Y
tree.corr <- function(EdgeP, Y, method, weight = F) {
  n <- ncol(EdgeP)
  
  if (method == "Dn") {
    phi <- apply(EdgeP, 1, function(P) independence::hoeffding.D.test(P, Y, precision = 1, collisions = F)$Dn)
  } else if (method == "Tn") {
    phi <- apply(EdgeP, 1, function(P) independence::tau.star.test(P, Y, precision = 1, collisions = F)$Tn)
  } else if (method == "Rn") {
    phi <- apply(EdgeP, 1, function(P) independence::hoeffding.refined.test(P, Y, precision = 1, collisions = F)$Rn)
  } else {
    stop("Please choose valid method.")
  }
  
  c(sum(weight * phi), max(phi))
}

tree.corr.cond <- function(EdgeP, Y, X, method, neighbor = min(10, length(Y)), weight = F) {
  n <- ncol(EdgeP)
  d <- nrow(EdgeP)
  distm <- as.matrix(stats::dist(X))
  phi <- rep(0, d)
  i.list <- apply(distm, 2, function(X) order(X, decreasing = F)[1:neighbor])
  Y1 <- matrix(Y[i.list], nrow = nrow(i.list))
  
  if (method == "Dn") {
    phi <- rowSums(apply(i.list, 2, function(x) apply(EdgeP[, x], 1, function(P) independence::hoeffding.D.test(P, Y[x], precision = 1, collisions = F)$Dn)))
  } else if (method == "Tn") {
    phi <- rowSums(apply(i.list, 2, function(x) apply(EdgeP[, x], 1, function(P) independence::tau.star.test(P, Y[x], precision = 1, collisions = F)$Tn)))
  } else if (method == "Rn") {
    phi <- rowSums(apply(i.list, 2, function(x) apply(EdgeP[, x], 1, function(P) independence::hoeffding.refined.test(P, Y[x], precision = 1, collisions = F)$Rn)))
  } else {
    stop("Please choose valid method.")
  }
  
  phi <- phi / n
  c(sum(weight * phi), max(phi))
}