#' Sparse community matrices
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
#' tree <- read.tree(text ="((t1:1,t2:1)N2:1,(t3:1,t4:1)N3:1)N1;")
#' com <- matrix(c(1,0,1,1,0,0,
#'                 1,0,0,1,1,0,
#'                 1,1,1,1,1,1,
#'                 0,0,1,1,0,1), 6, 4,
#'               dimnames=list(paste0("g",1:6), tree$tip.label))
#' com <- as.community(com)
#' com
#'
#' @rdname as.community
#' @export
as.community <- function (x, ...) {
  if (identical(class(x), "community"))
    return(x)
  UseMethod("as.community")
}


#' @rdname as.community
#' @export
as.community.matrix <- function(x){
  y <- vector("list", nrow(x))
  tmp <- seq_len(ncol(x))
  for(i in seq_len(nrow(x))){
    y[[i]] <- tmp[x[i,]>0]
  }
  names(y) <- rownames(x)
  attr(y, "labels") <- colnames(x)
  class(y) <- c("community", "splits")
  y
}


#' @rdname as.community
#' @export
read.community <- function(file, ...){
  d <- read.csv(file, ...)
#  M <- sparseMatrix(as.integer(d[,"grids"]), as.integer(d[,"species"]),
#                    dimnames = list(levels(d[,"grids"]), levels(d[,"species"])))
  M <- sparseMatrix(as.integer(d[,"species"]), as.integer(d[,"grids"]),
                    dimnames = list(levels(d[,"species"]), levels(d[,"grids"])))
  browser()
  res <- vector("list", ncol(M))
  for(i in seq_len(ncol(M))){
    res[[i]] <- M@i[(M@p[i]+1) : (M@p[i+1])] + 1
  }
  names(res) <- levels(d[,"grids"])
  attr(res, "labels") <- levels(d[,"species"])
  class(res) <- c("community", "splits")
  res
}


as.phylo_community <- function(tree, x){

}
