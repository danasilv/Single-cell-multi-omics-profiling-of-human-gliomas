library(phytools)
library(adephylo)

##################################################################### Transcriptional similarity ############################################################################
# tree.list is a list of all of the phylogenies (all tree replicates from all patient samples and plates) provided in this repository
# binmat is a matrix of gene module scores and contained within each element of tree.list
# used in Figure 5
#############################################################################################################################################################################

transcriptional.similarity <- function(tree, num_tests=10^4) {

  # gather pairwise expression and phylogenetic distances
  node_decay <- function(tree) {
        dat <- tree$binmat
        edif <- c(dist(dat))
        pdif <- c(distTips(tree, method="nNodes"))
        return(data.frame("expd"=edif, "node"=pdif))
  }
    x <- node_decay(tree)
    # measure correlation between transcriptional and lineage (node) distances
    obs <- cor(x$expd, x$node)
    ran <- rep(0, num_tests)
    # permutation test
    for(j in 1:num_tests) {
            tree0 <- tree
            tree0$binmat <- tree0$binmat[sample(1:length(tree0$tip.label), length(tree0$tip.label), replace = F),]
        xx <- node_decay(tree0)
        # correlation between transcriptional and lineage (node) distances in a randomly leaf-permuted tree
        ran[j] <- cor(xx$expd, xx$node)
    }
    # compute Z-score
    Z <- (obs - mean(ran))/sd(ran)
    return(Z)
}

# load phylogenies
tree.list <- readRDS("../treelist.RDS")

# output data frame of transcriptional similarity Z-scores
tsdf <- data.frame("bio.rep"=sapply(tree.list, '[[', "bio.rep"), "sample"=sapply(tree.list, '[[', "sample_name"),
                   "tree_rep"=sapply(tree.list, '[[', "tree_rep"), "class"=sapply(tree.list, '[[', "class"), "tran.sim"=sapply(tree.list, transcriptional.similarity))

fn <- paste0(paste("transcriptional_similarity_out", gsub(" ", "_", Sys.time() ) , sep = "_"), ".RDS")
saveRDS(tsdf, file=fn)
#
