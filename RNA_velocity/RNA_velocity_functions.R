library(Seurat)
library(SeuratWrappers)
library(velocyto.R)

######################### RNA Velocity Functions ##################################################################################################
# uses Seurat to process and integrate RNA velocity data, and SeuratWrappers and velocyto.R to run RNA velocity.
# ldat.idh and ldat.gbm are lists of the processed RNA velocity loom data each containing spliced, unspliced, ambiguous and spanning matrices.
# ldat.gbm, ldat.idh and the gene sets ("neftel_genes" and "tirosh_genes") need to be loaded from the repository before running code below
# executed using Seurat version 3.1.2, SeuratWrappers version 0.3.0 and velocyto.R version 0.6 in R 3.6.1
###################################################################################################################################################

# integrate data and run velocity in GBM
velocity.GBM <- function(plate.id, ldat) {

    s_obj <- list()
    # create a Seurat object for each plate.id provided using the supplied RNA velocity data (ldat)
    for (i in 1:length(plate.id))
    {
        keep <- grep(plate.id[i], colnames(ldat$unspliced))
        ldat_temp <- ldat
        ldat_temp$spliced <- ldat_temp$spliced[,keep]
        ldat_temp$unspliced <- ldat_temp$unspliced[,keep]
        ldat_temp$ambiguous <- ldat_temp$ambiguous[,keep]
        ldat_temp$spanning <- ldat_temp$spanning[,keep]
        s_obj[[i]] <- as.Seurat(ldat_temp)
    }

    # process each Seurat object and add gene-module scores (NPC-like, OPC-like, AC-like and MES-like)
    for (i in 1:length(s_obj)) {
        s_obj[[i]] <- NormalizeData(s_obj[[i]])
        s_obj[[i]] <- FindVariableFeatures(s_obj[[i]])
        s_obj[[i]] <- AddModuleScore(s_obj[[i]], neftel_genes, name=names(neftel_genes))
        s_obj[[i]] <- AddMetaData(s_obj[[i]], max.col(s_obj[[i]][[c("NPC1","OPC2", "AC3", "MES4")]]), col.name="cell_states")
    }

    # integrate Seurat objects
    s_anchors <- FindIntegrationAnchors(s_obj, dims = 1:15, k.filter = NA, k.score = 20)
    s_obj.int <- IntegrateData(s_anchors, dims = 1:15, preserve.order = TRUE)
    s_obj.int <- ScaleData(s_obj.int)
    s_obj.int <- RunPCA(s_obj.int,  features = c(neftel_genes$NPC, neftel_genes$OPC, neftel_genes$AC, neftel_genes$MES, cc.genes.updated.2019$s.genes, cc.genes.updated.2019$g2m.genes))

    # run RNA velocity
    s_obj.int <- RunVelocity(s_obj.int)

    return(s_obj.int)
}

# integrate data and run velcoity in IDH-MUT
velocity.IDH <- function(plate.id, ldat) {

    s_obj <- list()
    # create a Seurat object for each plate.id using the supplied RNA velocity data (ldat)
    for (i in 1:length(plate.id))
    {
        keep <- grep(plate.id[i], colnames(ldat$unspliced))
        ldat_temp <- ldat
        ldat_temp$spliced <- ldat_temp$spliced[,keep]
        ldat_temp$unspliced <- ldat_temp$unspliced[,keep]
        ldat_temp$ambiguous <- ldat_temp$ambiguous[,keep]
        ldat_temp$spanning <- ldat_temp$spanning[,keep]
        s_obj[[i]] <- as.Seurat(ldat_temp)
    }

    # process each Seurat object and add gene module scores (stem-like, AC-like and OC-like)
    for (i in 1:length(s_obj)) {
        s_obj[[i]] <- NormalizeData(s_obj[[i]])
        s_obj[[i]] <- FindVariableFeatures(s_obj[[i]])
        s_obj[[i]] <- AddModuleScore(s_obj[[i]], tirosh_genes, name=names(tirosh_genes))
        s_obj[[i]] <- AddMetaData(s_obj[[i]], max.col(s_obj[[i]][[c("stem1","AC2", "OC3")]]), col.name="cell_states")
    }

    # integrate Seurat objects
    s_anchors <- FindIntegrationAnchors(s_obj, dims = 1:15, k.filter = NA, k.score = 20)
    s_obj.int <- IntegrateData(s_anchors, dims = 1:15, preserve.order = TRUE)
    s_obj.int <- ScaleData(s_obj.int)
    s_obj.int <- RunPCA(s_obj.int,  features = c(tirosh_genes$stem, tirosh_genes$AC, tirosh_genes$OC, cc.genes.updated.2019$s.genes, cc.genes.updated.2019$g2m.genes))

    # run RNA velocity
    s_obj.int <- RunVelocity(s_obj.int)

    return(s_obj.int)
}

# simple Seurat pipeline to process GBM data and assign gene-module scores
simple_seu.gbm <- function(tpm) {
    x <- CreateSeuratObject(as.matrix(tpm))
    x <- NormalizeData(x)
    x <- ScaleData(x)
    x <- AddModuleScore(x, neftel_genes, name=names(neftel_genes), nbin = 30)
    x <- AddMetaData(x, max.col(x[[c("NPC1","OPC2", "AC3", "MES4")]]), col.name="neftel_states")
    return(x)
}

# estimate de-differentiation frequencies from RNA velocity in GBM
celltrans.gbm <- function(init.seu.obj, gbm.names) {

    # x[3] and x[4] refer to mature (AC and MES) module columns
    get.diff <- function(x) {
        z <- NA
        if((x[3] > 0 | x[4] > 0) ) {
            z <- 1
        }
        else {
            z <- 0
        }
        return(z)
    }

    # x[1] and x[2] refer to stem (NPC and OPC) module columns
    get.dediff <- function(x) {
        z <- NA
        if( (x[1] > 0 | x[2] > 0)) {
            z <- 1
        }
        else {
            z <- 0
        }
        return(z)
    }

    # simple_seu runs a basic Seurat pipeline on the RNA velocity projected expression patterns
    proj.seu.obj <- simple_seu.gbm(init.seu.obj@tools$RunVelocity$projected)

    # get initial and projected cell module scores
    init.mods <- init.seu.obj[[c("NPC1", "OPC2", "AC3", "MES4")]]
    proj.mods <- proj.seu.obj[[c("NPC1", "OPC2", "AC3", "MES4")]]

    # assign cell states
    init.seu.obj <- AddMetaData(init.seu.obj, as.factor(colnames(init.mods)[max.col(init.mods)]), col.name = "cell_states")
    proj.seu.obj <- AddMetaData(proj.seu.obj, as.factor(colnames(proj.mods)[max.col(proj.mods)]), col.name = "cell_states")

    # get difference between initial and projected modules scores
    mod.diffs <- proj.mods - init.mods

    out.mat <- matrix(NA, length(gbm.names), 1)
    rownames(out.mat) <- gbm.names
    colnames(out.mat) <- c("dediff.freq")
    init.seu.obj$sample <- NA

    for(j in 1:length(gbm.names)) {
        init.seu.obj$sample[grep(gbm.names[j], colnames(init.seu.obj))] <- gbm.names[j]
    }
    for(i in 1:length(gbm.names) ) {
        m1 <- as.matrix(mod.diffs[which(init.seu.obj$cell_states == "NPC1" & init.seu.obj$sample == gbm.names[i]),])
        m2 <- as.matrix(mod.diffs[which(init.seu.obj$cell_states == "OPC2" & init.seu.obj$sample == gbm.names[i]),])
        m3 <- as.matrix(mod.diffs[which(init.seu.obj$cell_states == "AC3"  & init.seu.obj$sample == gbm.names[i]),])
        m4 <- as.matrix(mod.diffs[which(init.seu.obj$cell_states == "MES4" & init.seu.obj$sample == gbm.names[i]),])

        # compute mature dediff frequency
        d2 <- (sum(apply(m3, 1, get.dediff)) + sum(apply(m4, 1, get.dediff)) )/(length(m3)+length(m4))

        out.mat[i,] <- d2
    }
    return(out.mat)
}
#

# simple Seurat pipeline to process IDH-MUT data and assign gene-module scores
simple_seu.idh <- function(tpm) {
    x <- CreateSeuratObject(as.matrix(tpm))
    x <- NormalizeData(x)
    x <- ScaleData(x)
    x <- AddModuleScore(x, tirosh_genes, name=names(tirosh_genes), nbin = 30)
    x <- AddMetaData(x, max.col(x[[c("stem1","AC2", "OC3")]]), col.name="tirosh_states")
    return(x)
}

# estimate de-differentiation frequencies from RNA velocity in IDH
celltrans.idh <- function(init.seu.obj, idh.names) {

    # x[2] and x[3] refer to mature (AC and OC) module columns
    get.diff <- function(x) {
        z <- NA
        if((x[2] > 0 | x[3] > 0) ) {
            z <- 1
        }
        else {
            z <- 0
        }
        return(z)
    }

    # x[1] refers to the stem module column
    get.dediff <- function(x) {
        z <- NA
        if( (x[1] > 0 )) {
            z <- 1
        }
        else {
            z <- 0
        }
        return(z)
    }

    # simple_seu runs a basic Seurat pipeline on the RNA Velocity projected expression patterns
    proj.seu.obj <- simple_seu.idh(init.seu.obj@tools$RunVelocity$projected)

    # get initial and projected cell module scores
    init.mods <- init.seu.obj[[c("stem1", "AC2", "OC3")]]
    proj.mods <- proj.seu.obj[[c("stem1", "AC2", "OC3")]]

    # assign cell states
    init.seu.obj <- AddMetaData(init.seu.obj, as.factor(colnames(init.mods)[max.col(init.mods)]), col.name = "cell_states")
    proj.seu.obj <- AddMetaData(proj.seu.obj, as.factor(colnames(proj.mods)[max.col(proj.mods)]), col.name = "cell_states")

    # get difference between initial and projected modules scores
    mod.diffs <- proj.mods - init.mods

    out.mat <- matrix(NA, length(idh.names), 1)
    rownames(out.mat) <- idh.names
    colnames(out.mat) <- c("dediff.freq")
    init.seu.obj$sample <- NA

    for(j in 1:length(idh.names)) {
        init.seu.obj$sample[grep(idh.names[j], colnames(init.seu.obj))] <- idh.names[j]
    }

    for(i in 1:length(idh.names)) {
        m1 <- as.matrix(mod.diffs[which(init.seu.obj$cell_states == "stem1" & init.seu.obj$sample == idh.names[i] ),])
        m2 <- as.matrix(mod.diffs[which(init.seu.obj$cell_states == "AC2" & init.seu.obj$sample == idh.names[i]),])
        m3 <- as.matrix(mod.diffs[which(init.seu.obj$cell_states == "OC3" & init.seu.obj$sample == idh.names[i]),])

        # compute mature dediff frequency
        d2 <- (sum(apply(m2, 1, get.dediff)) + sum(apply(m3, 1, get.dediff)) )/(length(m2)+length(m3))

        out.mat[i,] <- d2
    }
    return(out.mat)
}
#
