library(diversitree)
library(ape)
library(phytools)
library(GenSA)
library(bbmle)
library(adephylo)

################################# BiSSE Maximum Likelihood Estimation Pipeline ###############################################################
# tree.list is a list of all of the phylogenies (all tree replicates from all patient samples and plates) provided in this repository
# binmat is a matrix of gene module scores and contained within each element of tree.list
# trees within tree.list are annotated by class (GBM or IDH), patient.id (e.g. MGH105), sample_name (e.g. MGH105A), and tree_rep.
# used in Figure 5 and Extended Figure 10
###############################################################################################################################################

# load BiSSE functions
source("bisse_functions.R")

# load phylogenies
tree.list <- readRDS("../xcor/treelist.RDS")

# run BiSSE MLE on all phylogenies (note: this can take a while)
bisse.out <- list()
for(i in 1:length(tree.list)) {
  bisse.out[[i]] <- bisse.mle(tree.list[[i]])
  message("tree ", i, " complete")
}

# generate data frame of the most likely parameter estimates per phylogeny
bisse.out.df <- Reduce("rbind", lapply(bisse.out, get_best_coef))

fn <- paste0(paste("bisse_estimates", gsub(" ", "_", Sys.time() ) , sep = "_"), ".RDS")
saveRDS(bisse.out.df, file=fn)
#
