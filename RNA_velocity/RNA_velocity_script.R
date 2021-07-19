library(Seurat)
library(SeuratWrappers)
library(velocyto.R)


######################### RNA Velocity Pipeline ################################################################################################
# uses Seurat to process and integrate RNA velocity data, and SeuratWrappers and velocyto.R to run RNA velocity.
# ldat.idh and ldat.gbm are lists of the processed RNA velocity loom data each containing spliced, unspliced, ambiguous and spanning matrices.
# ldat.gbm, ldat.idh and the gene sets ("neftel_genes" and "tirosh_genes") need to be loaded from the repository before running code below
# executed using Seurat version 3.1.2, SeuratWrappers version 0.3.0 and velocyto.R version 0.6 in R 3.6.1
################################################################################################################################################

# load RNA velocity functions
source("RNA_velocity_functions.R")

# load data
ldat.gbm <- readRDS("ldat.GBM.RDS")
ldat.idh <- readRDS("ldat.IDH.RDS")

# load gene sets for GBM and IDH-MUT
neftel_genes <- readRDS("neftel_genes.RDS")
tirosh_genes <- readRDS("tirosh_genes.RDS")

# plate IDs for Seurat processing
gbm.nm.v <- c("MGH115", "MGH121_Plate_1", "MGH121_Plate_2", "MGH121_Plate_3", "MGH121_Plate_4", 
              "MGH122", "MGH124", "MGH105A", "MGH105B", "MGH105C", "MGH105D", "MGH211", "MGH129")
idh.nm.v <- c("MGH135", "MGH208.Pl1", "MGH208.Pl2", "MGH201", "MGH64", "MGH107", "MGH45", "MGH142.P.1", "MGH142.P.2")

# bio.rep IDs
gbm.nm.bio.rep <- c("MGH115", "MGH121", "MGH122", "MGH124", "MGH105A", "MGH105B", "MGH105C", "MGH105D", "MGH211", "MGH129")
idh.nm.bio.rep <- c("MGH135", "MGH208", "MGH201", "MGH64", "MGH107", "MGH45", "MGH142")

# for each tumor class, integrate data with Seurat and then run RNA velocity
vel.out.GBM <- velocity.GBM(gbm.nm.v, ldat.gbm)
vel.out.IDH <- velocity.IDH(idh.nm.v, ldat.idh)

# compute de-differentiation frequencies
dediff.vel.GBM <- celltrans.gbm(vel.out.GBM, gbm.nm.bio.rep)
dediff.vel.IDH <- celltrans.idh(vel.out.IDH, idh.nm.bio.rep)

vel_out <- rbind(dediff.vel.GBM, dediff.vel.IDH)
fn <- paste0(paste("velocity_out", gsub(" ", "_", Sys.time()),sep = "_"), ".RDS")
saveRDS(vel_out, file=fn)
#
