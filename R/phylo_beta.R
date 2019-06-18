# @importFrom betapart phylo.beta.multi
# @importFrom betapart phylo.beta.pair

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


phylo_community <- function(x, phy){
  el <- numeric(max(phy$edge))
  el[phy$edge[,2]] <- phy$edge.length
  x <- x[, phy$tip.label]
  if(is.matrix(x) | is(x, "sparseMatrix")){
    x <- as.splits(x)
  }
#  x <- phangorn::changeOrder(x, phy$tip.label)
  if(is.character(x) | is.numeric(x))  y <- list(phylo_com(x, phy))
  if(is.list(x)){
    y <- lapply(x, function(x, phy)phylo_com(x, phy), phy)
  }
  if(is.null(y)) return(NULL)
  attr(y, "edge.length") <- el
  attr(y, "labels") <- c(phy$tip.label, as.character(Ntip(phy)
                                                     + (1:Nnode(phy))))
  class(y) <- c("phylo_community", "splits")
  y
}


#' Phylogenetic beta diversity
#'
#' \code{beta_sparse_core} and \code{phylobeta_core} computes efficiently for
#' large community matrices and trees the nesseary quantitities used by the
#' betapart package to compute pairwise and multiple-site phylogenetic
#' dissimilarities.
#'
#' @aliases phylo_betapart_core beta_sparse_core
#' @param x an object of class Matrix or matrix
#' @param phy a phylogenetic tree (object of class phylo)
#'
#' @keywords cluster
#' @seealso read.community pd
#' @examples
#' library(ape)
#' tree <- read.tree(text ="((t1:1,t2:1)N2:1,(t3:1,t4:1)N3:1)N1;")
#' com <- matrix(c(1,0,1,1,0,0,
#'                 1,0,0,1,1,0,
#'                 1,1,1,1,1,1,
#'                 0,0,1,1,0,1), 6, 4,
#'               dimnames=list(paste0("g",1:6), tree$tip.label))
#' pbc <- phylobeta_core(com, tree)
#' library(betapart)
#' phylo.beta.multi(pbc)
#' phylo.beta.pair(pbc)
#'
#' @rdname phylobeta_core
#' @importFrom fastmatch fmatch
#' @export
phylobeta_core <- function(x, phy){
  x <- phylo_community(x, phy)
  l <- length(x)
  el <- attr(x, "edge.length")
  pd_tmp <- vapply(x, function(x, el)sum(el[x]), 0, el)
#  pd_tmp <- pd(x)

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
  res <- list(sumSi=sum(pd_tmp), St=sum(el), shared=SHARED,
              sum.not.shared = sum.not.shared,
              max.not.shared=max.not.shared, min.not.shared=min.not.shared)
  class(res) <- "phylo.betapart"
  res
}



phylobeta_core_old <- function(x){
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
  res <- list(sumSi=sum(pd_tmp), St=sum(el), shared=SHARED,
              sum.not.shared = sum.not.shared,
              max.not.shared=max.not.shared, min.not.shared=min.not.shared)
  class(res) <- "phylo.betapart"
  res
}

#' Match taxa and in phylogeny and community matrix
#'
#' match_phylo_comm compares taxa (species, labels, tips) present in a phylogeny
#' with a community matrix. Pruning, sorting and trying to add missing species
#' on genus level if possible to match in subsequent analysis.
#'
#' Based on the function of the same name in picante but allows sparse matrices
#' and with taxa addition.
#'
#' @param phy A phylogeny
#' @param comm A (sparse) community data matrix
#' @keywords cluster
#' @examples
#' example(phylo_builder)
#' pc <- match_phylo_comm(tree, comm)
#' @importFrom ape keep.tip
#' @export
match_phylo_comm <- function (phy, comm)
{
  if (!(is.data.frame(comm) | is.matrix(comm) | inherits(comm, "Matrix") )) {
    stop("Community data should be a data.frame or matrix with samples in rows
         and taxa in columns")
  }
  res <- vector("list", 2)
  commtaxa <- colnames(comm)
  if (is.null(commtaxa)) {
    stop("Community data set lacks taxa (column) names, these are required to
         match phylogeny and community data")
  }
  phytaxa <- phy$tip.label
  index <- intersect(commtaxa, phytaxa)
  res$comm <- comm[, index]
  res$phy <- keep.tip(phy, index)
  return(res)
}


#' @rdname phylobeta_core
#' @importFrom Matrix Matrix tcrossprod colSums
#' @export
## non_phylo version
beta_sparse_core <- function (x) {
  if (!inherits(x, "Matrix")) x <- Matrix(x)

  shared <- as.matrix( tcrossprod(x) )# %*% t(x)
  not.shared <- abs(sweep(shared, 2, diag(shared)))
  sumSi <- sum(diag(shared))
  St <- sum(colSums(x) > 0)
  a <- sumSi - St
  sum.not.shared <- not.shared + t(not.shared)
  max.not.shared <- pmax(not.shared, t(not.shared))
  min.not.shared <- pmin(not.shared, t(not.shared))
  computations <- list(data = x, sumSi = sumSi, St = St, a = a, shared = shared,
                       not.shared = not.shared, sum.not.shared = sum.not.shared,
                       max.not.shared = max.not.shared,
                       min.not.shared = min.not.shared)
  class(computations) <- "betapart"
  return(computations)
}
