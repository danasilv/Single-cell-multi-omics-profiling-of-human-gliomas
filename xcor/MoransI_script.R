library(phytools)
library(adephylo)
library(dplyr)

################################################################## Moran's I Permutation Test ##############################################################################################
# used in Figure 4 and Extended Figure 10
############################################################################################################################################################################################

# load xcor functions
source("xcor_functions.R")

# load phylogenies
tree.list <- readRDS("treelist.RDS")

# GBM tree replicate indices in tree.list
id1 <- 1:110

# IDH-MUT tree replicate indices in tree.list
id2 <- 111:160

# get Moran's I, p-value, and Z-score from a 10^6 leaf-permutation test for GBM
M1 <- data.frame(t(sapply(tree.list[id1], function(x) {m <- get.moran(x$binmat, proxTips(x, method="nNodes", normalize = "none"), nrepet = 10^6); return(c(m$obs, m$pvalue, m$Z)) })),
                "bio.rep"=sapply(tree.list[id1], '[[', "bio.rep"),
                "sample_name"=sapply(tree.list[id1], '[[', "sample_name"),
                "tree_rep"=sapply(tree.list[id1], '[[', "tree_rep"), 
                "class"=sapply(tree.list[id1], '[[', "class"))

colnames(M1)[1:12] <- c("GBM.NPC.moran.I", "GBM.OPC.moran.I", "GBM.AC.moran.I", "GBM.MES.moran.I", "GBM.NPC.moran.P", "GBM.OPC.moran.P",
                        "GBM.AC.moran.P", "GBM.MES.moran.P", "GBM.NPC.moran.Z", "GBM.OPC.moran.Z", "GBM.AC.moran.Z", "GBM.MES.moran.Z")

# get Moran's I, p-value, and Z-score from a 10^6 leaf-permutation test for IDH-MUT
M2 <- data.frame(t(sapply(tree.list[id2], function(x) {m <- get.moran(x$binmat, proxTips(x, method="nNodes", normalize = "none"), nrepet = 10^6); return(c(m$obs, m$pvalue, m$Z)) })),
                "bio.rep"=sapply(tree.list[id2], '[[', "bio.rep"),
                "sample_name"=sapply(tree.list[id2], '[[', "sample_name"),
                "tree_rep"=sapply(tree.list[id2], '[[', "tree_rep"), 
                "class"=sapply(tree.list[id2], '[[', "class"))

colnames(M2)[1:9] <- c("IDH.Stem.moran.I", "IDH.AC.moran.I", "IDH.OC.moran.I", "IDH.Stem.moran.P", "IDH.AC.moran.P", "IDH.OC.moran.P",
                       "IDH.Stem.moran.Z", "IDH.AC.moran.Z", "IDH.OC.moran.Z")


ti <- Sys.time()
fn1 <- paste0(paste("GBM_MoransI", gsub(" ", "_", ti) , sep = "_"), ".RDS")
fn2 <- paste0(paste("IDH_MoransI", gsub(" ", "_", ti) , sep = "_"), ".RDS")
saveRDS(M1, file=fn1)
saveRDS(M2, file=fn2)
#
