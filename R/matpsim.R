#' A function to calculate phylogenetic beta diversity.
#'
#' A function to calculate phylogenetic beta diversity.
#'
#' Using Simpson's phylogenetic beta diversity
#'
#' @param phyl Phylogenetic tree
#' @param com Distribution data, can be point records or range maps, represented as species × site matrix
#' @param clust Number of computing cores to use. Defaults to 2
#' @keywords cluster
#' @export
#' @examples
#' tree <- read.tree(text ="((t1:1,t2:1)N2:1,(t3:1,t4:1)N3:1)N1;")
#' com <- matrix(c(1,0,1,1,0,0,
#'                 1,0,0,1,1,0,
#'                 1,1,1,1,1,1,
#'                 0,0,1,1,0,1), 6, 4,
#'                 dimnames=list(paste0("g",1:6), tree$tip.label))
#' matpsim(tree, com, 1L)



matpsim <- function(phyl, com, clust = 2) # make sure nodes are labelled and that com and phyl species match. Returns as dist object with pairwise phylobetasim distances
{

#  require(phylobase)  # detailed information for all phylo branches
  new <- phylo4(phyl)
  dat <- as.data.frame(print(new))
  allbr <- dat$edge.length
  names(allbr) <- getEdge(new)
  spp <- colnames(com)

 #   require(foreach)



  # create a list of phylo branches for each species

  brs <-  foreach(i = spp, .packages = "phylobase") %do% #this loop makes a list of branches for each species
  {
    print(which(spp == i)/length(spp))
    print(date())
    brsp <- vector()
    br   <- as.numeric(rownames(dat[dat$label == i,]))
    repeat{
      brsn <- getEdge(new,br)
      brsl <- dat[names(brsn),"edge.length"]
      names(brsl) <- brsn
      brsp <- c(brsp, brsl)
      br   <- dat[br,3]
      if(br == 0) {break}
    }

    brsp
  }
  names(brs) <- spp




  print("brs")

  # create a species by phy branch matrix

  spp_br <- matrix(0,nrow = length(spp), ncol = length(allbr))
  rownames(spp_br) <- spp
  colnames(spp_br) <- names(allbr)

  for(i in spp)
  {
    spp_br[i,names(brs[[i]])] <- brs[[i]]
  }

  spp_br <- spp_br[,-(ncol(com)+1)] # removes root
  spp_br <- spp_br[,!colSums(spp_br) %in% c(0,nrow(spp_br))]  # here take out all common branches instead

  print("spp_br")

  spp_br <<- spp_br  ### this send the species by phy branch matrix to the workspace



  # function to give the lengths of phy branches within cell i

  cellbr <- function(i,spp_br, com)
  {
    i_spp <- rownames(spp_br)[com[i,]>0]
    if(length(i_spp) > 1)
    {i_br  <- as.numeric(apply(spp_br[i_spp,],2,max))}
    else  {i_br  <- as.numeric(spp_br[i_spp,]) }
    names(i_br) <- colnames(spp_br)
    return(i_br)
  }


  #require(foreach)
  #require(doSNOW)



  # create cell by phy branch matrix, use branch lengths rather than presence/absense

  tcellbr <- foreach(j = rownames(com), .combine = "rbind") %do% {cellbr(j,spp_br,com)}

  tcellbr <- tcellbr[,!is.na(colSums(tcellbr))] # new line


  print("cell_br")

  rownames(tcellbr) <- rownames(com)

  tcellbr <<- tcellbr ## send the cell by phy branch matrix to the workspace

  # function to calculate phylobsim between cell_a and cell_b

  nmatsim <- function(cell_b) # samp = grid cell of interest
  {
    a_br  <- tcellbr[cell_a,]
    b_br <- tcellbr[cell_b,]
    s_br <- rbind(a_br,b_br)
    s_br <- as.matrix(s_br[,colSums(s_br) >0])
    pa_br <- s_br > 0
    both <- s_br[1,colSums(pa_br > 0)==2]
    ubr <- s_br[,colSums(pa_br > 0)==1]
    a_ubr <- as.matrix(ubr)[1,]
    b_ubr <- as.matrix(ubr)[2,]
    psim <- 1 - (sum(both)/(min(sum(a_ubr),sum(b_ubr))+sum(both)))
    return(psim)
  }

  # calculate full cell by cell phylobsim matrix





  psim <- foreach(j = rownames(tcellbr)) %do% {print(j);cell_a <- j; unlist(lapply(rownames(tcellbr)[1:(which(rownames(tcellbr) == j))], nmatsim))}



  psim <- do.call("rbind", psim)
  rownames(psim) <- rownames(tcellbr)
  colnames(psim) <- rownames(tcellbr)
  psim[upper.tri(psim)] <- psim[lower.tri(psim)]
  psim <- as.dist(psim)
  return(psim)

}
