#' Create a subtree with largest overlap from a species list.
#'
#' Create a subtree with largest overlap from a species list.
#'
#' If species is not in the tip label species is added to the genus or family
#' level when possible
#'
#' @param splist Species list
#' @param tree a phylogenetic tree (object of class phylo)
#' @param nodes node list
#' @param output.splist return species list
#' @keywords cluster
#' @export
phylo_maker <- function (splist, tree, nodes, output.splist = TRUE)
{
    splist[sapply(splist, is.factor)] <-
            lapply(splist[sapply(splist, is.factor)], as.character)
    nodes[sapply(nodes, is.factor)] <-
             lapply(nodes[sapply(nodes, is.factor)], as.character)
    if (any(duplicated(splist$species))){
        warning("Duplicated species detected and removed.")
        splist <- splist[!duplicated(splist$species), ]
    }
    Ori <- splist
    splist$species <- gsub(" ", "_", splist$species)
    splist$species <- gsub("(^[[:alpha:]])", "\\U\\1", splist$species,
                           perl = TRUE)
    splist$genus <- gsub("(^[[:alpha:]])", "\\U\\1", splist$genus,
                         perl = TRUE)
    splist$family <- gsub("(^[[:alpha:]])", "\\U\\1", splist$family,
                          perl = TRUE)
    renameNode <- data.frame(node.label = paste("N", 1:length(tree$node.label),
                           sep = ""), original.node.label = tree$node.label,
                           stringsAsFactors = FALSE)
    tree$node.label <- paste("N", 1:length(tree$node.label),
                             sep = "")
    add.tip <- splist[which(is.na(match(splist$species, tree$tip.label))), ]
    status <- rep("match(prune)", dim(splist)[1])
    status[which(is.na(match(splist$species, tree$tip.label)))] <- "match(add)"
    if (dim(add.tip)[1] == 0)
        stop("Format mismatch OR There is no species needs to be added.")
    add.tip$sort <- ""
    add.tip$sort[which(!is.na(match(add.tip$genus,
                                nodes[nodes$level == "G", ]$genus)))] <- "G"
    add.tip$sort[which(is.na(match(add.tip$genus, nodes[nodes$level == "G", ]$genus))
          & !is.na(match(add.tip$family, nodes[nodes$level == "F", ]$family)))] <- "F1"
    add.tip$sort[add.tip$sort == "F1"][duplicated(add.tip[add.tip$sort ==
                                                        "F1", ]$genus)] <- "F2"

    a <- which(add.tip$sort == "")
    b <- as.character(add.tip$species[a])
    if (length(a) > 0)
        print(paste("Note:", length(a), "taxa unmatched."))
    status[match(b, splist$species)] <- "unmatch"
    Ori$status <- status
#    if ("S1" %in% scenarios) {

        T1 <- tree
        nodeG <- nodes[nodes$level == "G", ]
        nodeF <- nodes[nodes$level == "F", ]
        data <- add.tip[add.tip$sort == "F1", ]

        if (dim(data)[1] > 0) {
            n <- match(data$family, nodeF$family)
            num <- match(nodeF$node.label[n], T1$node.label) + length(T1$tip.label)
            T1 <- add.tips(T1, tips=data$species, num)
        }
        data <- add.tip[add.tip$sort == "F2", ]
        if (dim(data)[1] > 0) {

            nn <- lapply(paste(data$genus, "_", sep = ""), grep, T1$tip.label, fixed=TRUE)
            ll <- lengths(nn)
            num <- numeric(length(nn))
            if(any(ll==1)){
#                browser()
                tmp_num <- unlist(nn[ll==1])
                num <- T1$edge[match(tmp_num, T1$edge[,2]), 1]
                T1 <- add.tips(T1, tips=data$species[ll==1], where = num)
            }
            if(any(ll > 1)){
                ind <- which(ll>1)
                for (i in ind) {
                    n <- grep(paste(data$genus[i], "_", sep = ""),
                          T1$tip.label)
                    num <- getMRCA(T1, T1$tip.label[ nn[[i]] ])
                    T1 <- add.tips(T1, data$species[i], where = num)
                }
            }
        }
        data <- add.tip[add.tip$sort == "G", ]
        if (dim(data)[1] > 0) {
            n <- match(data$genus, nodeG$genus)
            num <- match(nodeG$node.label[n], T1$node.label) + length(T1$tip.label)
            T1 <- add.tips(T1, tips=data$species, num)
        }
        #        T1$edge.length <- round(T1$edge.length, 5)
        toDrop <- setdiff(1:length(T1$tip.label),
                          which(!is.na(match(T1$tip.label, splist$species))))
        T1 <- drop.tip(T1, tip = toDrop)
        Re <- which(!is.na(match(T1$node.label, renameNode$node.label)))
        noRe <- which(is.na(match(T1$node.label, renameNode$node.label)))
        T1$node.label[Re] <- renameNode$original.node.label[match(T1$node.label,
                                                    renameNode$node.label)[Re]]
        T1$node.label[noRe] <- ""

    if (output.splist == FALSE)
        splist <- NULL
    phylo <- list(Scenario.1 = T1, Species.list = Ori)
    phylo[sapply(phylo, is.null)] <- NULL
    return(phylo)
}


