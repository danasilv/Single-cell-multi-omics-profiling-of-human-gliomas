args = commandArgs(TRUE)

library(Homo.sapiens)

## List BAMs in Directory and Path to annotation file (.GTF)
expr_files <- list.files(path = snakemake@config[["path"]], pattern = "*.RSEM.genes.results$")
expr_matrices <- lapply(expr_files, read.table, sep = "\t", row.names = 1, header = TRUE)
expr_matrices <- lapply(expr_matrices, subset, select = "TPM")
expr_matrix <- Reduce(function(df1, df2) cbind(df1, df2), expr_matrices)
cell_names <- sapply(expr_files, function(x) strsplit(x, ".", fixed = TRUE)[[1]][1])
colnames(expr_matrix) <- cell_names

## Extract the gene SYMBOL and ENSEMBL ID given ENTREZID accession
info <- select( Homo.sapiens , keys = rownames(expr_matrix), columns = c("SYMBOL", "ALIAS"), keytype = "ENSEMBL")

## Remove Duplicated Gene SYMBOLs
expr_matrix <- expr_matrix[!duplicated(info$SYMBOL[match(rownames(expr_matrix), info$ENSEMBL) ]), ]

## Remove NA's
expr_matrix <- expr_matrix[!is.na(info$SYMBOL[match(rownames(expr_matrix), info$ENSEMBL) ]), ]

## Replace the rownames of expr_matrix (which are ENTREZID ACCs) with the gene's SYMBOL
rownames(expr_matrix) <- info$SYMBOL[match(rownames(expr_matrix), info$ENSEMBL)]

## Write out the expression matrix
write.table(expr_matrix, snakemake@output[["matrix"]], col.names = TRUE, row.names = TRUE, quote = FALSE, sep = "\t")
