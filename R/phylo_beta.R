#
# occurence --> phylogenetic occurence --> phylo beta
#           --> beta
#  beta  == beta core (a,b,c)
#  Matrix package, sparse binary data
#
#  betapart indices Jaccard, Simpson, etc.
#



phylo_com_C <- function(phy, tip){
  if (!inherits(phy, "phylo"))
    stop("object \"phy\" is not of class \"phylo\"")
  Ntips <- length(phy$tip.label)
  Nnodes <- phy$Nnode
  rootnd <- Ntips +1L #getRoot(phy)
  if (is.character(tip)) tip <- match(tip, c(phy$tip.label, phy$node.label))
  tip <- as.integer(tip)
  pvec <- integer(max(phy$edge))
  pvec[phy$edge[, 2]] <- phy$edge[, 1]
  res <- phylo_com_Rcpp(tip, pvec, Ntips, Nnodes, rootnd)
  res
}

phylo_com <- function(phy, tip){
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


# phylogenetic diversity
pd <- function(x, tree){
  el <- numeric(max(tree$edge))
  el[tree$edge[,2]] <- tree$edge.length
#  fun <- function(x, tree){
#    el <- numeric(max(tree$edge))
#    el[tree$edge[,2]] <- tree$edge.length
#    el[x]
#  }
  if(is.list(x)) lapply(x, function(x, el)el[x], el)
  el[x] #fun(x, tree)
}



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



cppFunction("IntegerVector phylo_com_Rcpp(IntegerVector tip, IntegerVector pvec,
             int ntip, int nnode, int rootnd) {
            IntegerVector result = clone(tip);
            int l = tip.size();
            result.push_back(ntip); //?
            IntegerVector done_v(ntip + nnode + 1); // ntip + nnode + 1
            done_v[rootnd] = 1;
            int nd=0, done=0;
            for(int i = 0; i<l; i++){
              nd = pvec[result[i]-1];
              done = done_v[nd];
              while(done < 1){
                done_v[nd] = 1;
                result.push_back(nd);
                nd = pvec[nd - 1];
                done = done_v[nd];
              }
            }
            std::sort(result.begin(), result.end());
            return result;
            }")



