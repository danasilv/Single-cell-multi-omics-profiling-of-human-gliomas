library(phytools)
library(adephylo)
library(ggplot2)
library(ggpubr)
library(cowplot)
library(dplyr)
library(reshape2)
library(Seurat)

############################################################## Cross-Correlation Functions #####################################################################
# contains functions to compute auto- and cross-correlation (univariate and bivariate Moran's I)
# applied in analyses shown in Figures 4 & 5 and Extended Figure 10
################################################################################################################################################################

# Moran's I Z-score calculated from moments
moran.analytical <- function(data, weight.matrix) {
    # data format: one row per cell and one column per measurement (e.g. one column per gene expression score)
    # Moment calculations from Czaplewski and Reich 1993.

    # from phylogenetic pairwise distances
    W <- weight.matrix
    # number of leaves
    n <- nrow(data)

    # center data
    d0 <- apply(data, 2, function(xx) xx - mean(xx))

    d1 <- t(d0)%*%d0/n
    d1.2 <- diag(d1)%*%t(diag(d1))
    d2 <- t(d0^2)%*%(d0^2)/n

    # below uses similar notation as in Czaplewski and Reich 1993
    x <- rep(1, n)
    w <- as.numeric(t(x)%*%W%*%x)
    S3 <- t(x)%*%(W * t(W))%*%x
    S4 <- t(x)%*%(W * W)%*%x
    S5 <- t(x)%*%W%*%W%*%x
    S6 <- t(W%*%x)%*%(W%*%x) + t(t(W)%*%x)%*%(t(W)%*%x)
    S1 <- S4 + S3
    S2 <- 2*S5 + S6

    S1 <- as.numeric(S1)
    S2 <- as.numeric(S2)
    S3 <- as.numeric(S3)
    S4 <- as.numeric(S4)
    S5 <- as.numeric(S5)
    S6 <- as.numeric(S6)

    # Moran's I 
    I <- (t(d0)%*%W%*%d0)/(w*sqrt(d1.2))

    # Expected value of Moran's I
    E.I2 <- -cor(data)/(n-1)

    # Variance of Moran's I
    V00 <- (  ((d1^2 * n)/ (d1.2) ) * (2*(w^2 - S2 + S1) + (2*S3 - 2*S5)*(n-3) + S3*(n-2)*(n-3) )  +  (( -d2 )/( d1.2 ))*(6*(w^2 - S2 + S1) + (4*S1-2*S2)*(n-3) + S1*(n-2)*(n-3)  ) 
              + n*((w^2 - S2 + S1) + (2*S4 -S6)*(n-3) + S4*(n-2)*(n-3) )  ) / ((n-1)*(n-2)*(n-3)*(w^2) ) -  (  ( ( d1^2 )/( d1.2 ) ) * ( 1/((n-1)^2)  )  )

    # Z-score calculation: (Measured I - Expected Value)/sqrt(Variance)
    Z0 <- (I-E.I2)/sqrt(V00)

    return(list("Morans.I"=I, "Expected.I"=E.I2,  "Var.I"=V00,"Z.score"=Z0))
}
#

# function to calculate (analytical) Moran's I taking a phylogeny as input
# feature should be set to data type (e.g. binmat is matrix of gene module scores)
# meth is the type of tree distance to use; e.g. "nNodes" refers to node distance
msig.tree <- function(tree, feature="binmat", meth="nNodes") {
  # get weight matrix
  W <- adephylo::proxTips(tree, method=meth, normalize="none")
  moran.analytical(tree[[feature]], W)
}
#

# computes mean Moran's I Z-score matrix given a vector of "sample" names
get_xcor_matrix <- function(sample) {

    l <- list()
    for(i in 1:length(sample)) {
      # searches tree.list for phylogeny replicates corresponding to selected sample
      r <- grep(sample[i], sapply(tree.list, '[[', "tree_rep"))

      # computes mean Z-score for each sample
      l[[i]] <- apply(simplify2array(lapply(tree.list[r], function(x) msig.tree(x)$Z)), 1:2, mean)

    }
    # computes mean Z-score across samples
    x <- apply(simplify2array(l), 1:2, mean)

    return(x)
}

# get.moran is a modified version of the function adephylo::abouheif.moran to output Z-scores
# used for leaf-permutation tests
get.moran <- function (x, W = NULL, method = c("oriAbouheif", "patristic", 
                                                "nNodes", "Abouheif", "sumDD"), f = function(x) {
                                                    1/x
                                                }, nrepet = 999, alter = c("greater", "less", "two-sided")) 
{
    alter <- match.arg(alter)
    method <- match.arg(method)
    if (!is.null(W)) {
        if (any(W < 0)) 
            stop("negative terms found in 'W'")
        if (nrow(W) != ncol(W)) 
            stop("'W' is not squared")
        W <- as.matrix(W)
    }
    else {
        if (!inherits(x, "phylo4d")) 
            stop("if W is not provided, x has to be a phylo4d object")
        if (is.character(chk <- checkPhylo4(x))) 
            stop("bad phylo4d object: ", chk)
        W <- proxTips(x, method = method, f = f, normalize = "row", 
                      symmetric = TRUE)
    }
    nobs <- ncol(W)
    W <- (W + t(W))/2
    if (inherits(x, "phylo4d")) {
        if (is.character(chk <- checkPhylo4(x))) 
            stop("bad phylo4d object: ", chk)
        x <- tdata(x, type = "tip")
    }
    x <- data.frame(x)
    test.names <- names(x)
    x <- data.matrix(x)
    if (nrow(x) != nobs) 
        stop("non convenient dimension")
    nvar <- ncol(x)
    res <- .C("gearymoran", param = as.integer(c(nobs, nvar, 
                                                 nrepet)), data = as.double(x), W = as.double(W), obs = double(nvar), 
              result = double(nrepet * nvar), obstot = double(1), restot = double(nrepet), 
              PACKAGE = "adephylo")
    res <- as.krandtest(obs = res$obs, sim = matrix(res$result, 
                                                    ncol = nvar, byrow = TRUE), names = test.names, alter = alter, output = "full")
    Z <- (res$obs - mean(res$sim))/sd(res$sim) # added Z-score calculation
    res$Z <- Z
    return(res)
}
#

# get pairwise cross-correlations and correlations from a list of phylogenies
# trees within tree_list should have the same number of genes and be in the same order
get_xcor_by_cor_df <- function(tree_list, genes) {

    # compute analytical cross-correlation Z-scores for each phylogeny
    xl <- list()
    for(i in 1:length(tree_list)) {
        xl[[i]] <- msig.tree(tree_list[[i]], feature = genes)$Z.score
    }

    # compute mean Z-scores across phylogenies
    xx <- apply(simplify2array(xl), 1:2, mean)

    # omit autocorrelations and construct data frame
    xx[upper.tri(xx, diag = T)] <- NA
    df1 <- reshape2::melt(xx, na.rm = T)

    # compute correlations
    rr <- cor(tree_list[[1]][[genes]])

    # omit self-correlations and construct data frame
    rr[upper.tri(rr, diag = T)] <- NA
    df2 <- reshape2::melt(rr, na.rm = T)

    colnames(df1) <- c("gene1", "gene2", "xcor")
    colnames(df2)[3] <- "cor"
    df1$cor <- df2$cor

    df1$cc.ident <- "none"
    df1$stem.ident <- "none"

    # get set of cycling genes
    cyc00 <- unique(c(cc.genes.updated.2019$s.genes, cc.genes.updated.2019$g2m.genes))

    # get GBM stem-like genes
    gbm_stem_genes <- unlist(neftel_genes[c("NPC", "OPC" )], use.names = F)

    # annotate gene pairs that are either both in the cycling or stem-like gene set
    df1[ which( (df1$gene1 %in% cyc00) & (df1$gene2 %in% cyc00 ) ), "cc.ident"] <- "cycling"
    df1[ which( (df1$gene1 %in% gbm_stem_genes ) & (df1$gene2 %in% gbm_stem_genes ) ), "stem.ident"] <- "stem"
    df1$cc.ident <- factor(df1$cc.ident, levels = c("none", "cycling"), ordered = TRUE)
    df1$stem.ident <- factor(df1$stem.ident, levels = c("none", "stem"), ordered = TRUE)

    return(df1)
}
#
