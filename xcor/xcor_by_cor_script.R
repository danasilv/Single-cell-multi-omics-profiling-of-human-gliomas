library(phytools)
library(adephylo)
library(ggplot2)
library(ggpubr)
library(cowplot)
library(dplyr)
library(reshape2)
library(Seurat)

################## Cross-Correlation by Correlation Script ######################
# used in Figure 4 and Extended Figure 10
#################################################################################

# load xcor functions
source("xcor_functions.R")

# load GBM gene sets
neftel_genes <- readRDS("../RNA_velocity/neftel_genes.RDS")

# load MGH115 and MGH122 phylogeny lists
mgh115.list <- readRDS("mgh115-list.RDS")
mgh122.list <- readRDS("mgh122-list.RDS")

# compute xcor and cor
xcorXcor115 <- get_xcor_by_cor_df(mgh115.list, "tpm")
xcorXcor122 <- get_xcor_by_cor_df(mgh122.list, "tpm")

ti <- Sys.time()
fn1 <- paste0(paste("xcorXcor115", gsub(" ", "_", ti) , sep = "_"), ".RDS")
fn2 <- paste0(paste("xcorXcor122", gsub(" ", "_", ti) , sep = "_"), ".RDS")
saveRDS(xcorXcor115, file=fn1)
saveRDS(xcorXcor122, file=fn2)
#
