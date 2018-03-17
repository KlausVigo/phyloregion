#
# occurence --> phylogenetic occurence --> phylo beta
#           --> beta
#  beta  == beta core (a,b,c)
#  Matrix package, sparse binary data
#
#  betapart indices Jaccard, Simpson, etc.
#

phylo_com <- function(tip, phy){
    if (!inherits(phy, "phylo"))
        stop("object \"phy\" is not of class \"phylo\"")
    Ntips <- length(phy$tip.label)
    Nnodes <- phy$Nnode
    done_v <- logical(Ntips + Nnodes)
    rootnd <- Ntips +1L #getRoot(phy)
    if (is.character(tip)) tip <- match(tip, c(phy$tip.label, phy$node.label))
    tip <- as.integer(tip)
    res <- pvec <- integer(max(phy$edge))
    pvec[phy$edge[, 2]] <- phy$edge[, 1]

    res[seq_along(tip)] <- tip
    l <- length(tip) + 1L
    res[l] <- rootnd
    done_v[rootnd] <- TRUE
    for(k in tip) {
        nd <- pvec[k]
        done <- done_v[nd]
        while(!done){
            done_v[nd] <- TRUE
            l <- l+1
            res[l] <- nd
            nd <- pvec[nd]
            done <- done_v[nd]
        }
    }
    sort(res[1:l])
}

#' Phylogenetic community objects
#'
#' read.community reads in file containing occurence data and returns a sparse
#' matrix.
#'
#' If species is not in the tip label species is added to the genus or family
#' level when possible
#'
#' @aliases read.community as.community as.community.matrix print.community
#' @param splist Species list
#' @param tree a phylogenetic tree (object of class phylo)
#' @param nodes node list
#' @param output.splist return species list
#' @keywords cluster
#' @examples
#' example(as.community)
#' pc <- phylo_community(com, tree)
#' pd(pc)
#'
#' @rdname phylo_community
#' @export
phylo_community <- function(x, phy){
  el <- numeric(max(phy$edge))
  el[phy$edge[,2]] <- phy$edge.length
  if(is.character(x) | is.numeric(x))  y <- list(phylo_com(x, phy))
  if(is.list(x)){
    y <- lapply(x, function(x, phy)phylo_com(x, phy), phy)
  }
  if(is.null(y)) return(NULL)
  attr(y, "edge.length") <- el
  attr(y, "labels") <- c(phy$tip.label, as.character(Ntip(phy) + (1:Nnode(phy))))
  class(y) <- c("phylo_community", "splits")
  y
}


#' @rdname phylo_community
#' @export
pd <- function(x, tree=NULL){
  if(!is.null(tree)){
    el <- numeric(max(tree$edge))
    el[tree$edge[,2]] <- tree$edge.length
  }
  else el <- attr(x, "edge.length")
  if(is.list(x)) res <- sapply(x, function(x, el)sum(el[x]), el)
  else res <- sum(el[x]) #fun(x, tree)
  res
}


#' @rdname phylo_community
#' @export
phylo_betapart_core <- function(x){
  l <- length(x)
  pd_tmp <- pd(x)
  el <- attr(x, "edge.length")
  Labels <- names(x)
  class(x) <- NULL
  SHARED <- vector("numeric", l*(l-1)/2)
  B <- vector("numeric", l*(l-1)/2)
  C <- vector("numeric", l*(l-1)/2)

  k <- 1
  for(i in 1:(l-1)){
    xi <- x[[i]]
    for(j in (i+1):l){
      #   sum(el[fast_intersect(x[[i]], x[[j]])])
      SHARED[k] <- sum(el[xi[fmatch(x[[j]], xi, 0L)] ])
      B[k] <- pd_tmp[i] - SHARED[k]
      C[k] <- pd_tmp[j] - SHARED[k]
      k <- k+1
    }
  }

  sum.not.shared <- B+C
  max.not.shared <- pmax(B,C)
  min.not.shared <- pmin(B,C)

  at <- structure(list(Labels=Labels, Size = l, class = "dist", Diag = FALSE,
                       Upper = FALSE), .Names = c("Labels", "Size", "class", "Diag", "Upper"))
  attributes(SHARED) <- at
  attributes(sum.not.shared) <- at
  attributes(max.not.shared) <- at
  attributes(min.not.shared) <- at
  res <- list(sumSi=sum(pd_tmp), St=sum(el), shared=SHARED, sum.not.shared = sum.not.shared,
              max.not.shared=max.not.shared, min.not.shared=min.not.shared)
  class(res) <- "phylo.betapart"
  res
}

# based on picante needs improvement
match.phylo.comm <- function (phy, comm, trace=1)
{
  if (!(is.data.frame(comm) | is.matrix(comm) | inherits(comm, "Matrix") )) {
    stop("Community data should be a data.frame or matrix with samples in rows and taxa in columns")
  }
  res <- list()
  phytaxa <- phy$tip.label
  commtaxa <- colnames(comm)
  if (is.null(commtaxa)) {
    stop("Community data set lacks taxa (column) names, these are required to
         match phylogeny and community data")
  }
  if (!all(commtaxa %in% phytaxa)) {
    if(trace){
      print("Dropping taxa from the community because they are not present in
            the phylogeny:")
      print(setdiff(commtaxa, phytaxa))
    }
    comm <- comm[, intersect(commtaxa, phytaxa)]
    commtaxa <- colnames(comm)
  }
  if (any(!(phytaxa %in% commtaxa))) {
    if(trace){
      print("Dropping tips from the tree because they are not present in
                    the community data:")
      print(setdiff(phytaxa, commtaxa))
    }
    res$phy <- prune.sample(comm, phy)
  }
  else {
    res$phy <- phy
  }
  res$comm <- comm[, res$phy$tip.label]
  return(res)
}

