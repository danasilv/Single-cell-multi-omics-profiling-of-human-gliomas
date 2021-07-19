library(Seurat)
library(dplyr)
library(Homo.sapiens)
library(mygene)
library(OneR)

## Load in Tirosh 2016 modules (IDH)
tirosh = read.table("Tirosh2016_metamodules.txt", sep = '\t', header = T, na.strings = "")
names(tirosh) <- c("G1.S", "G2.M", "stemness", "AC_PCA", "AC_PCA_MICE", "OC_PCA", "OC_PCA_MICE")
tirosh$G1.S = gsub('\\s+', '', tirosh$G1.S)
tirosh$G2.M = gsub('\\s+', '', tirosh$G2.M)
tirosh$stemness = gsub('\\s+', '', tirosh$stemness)
tirosh$AC_PCA = gsub('\\s+', '', tirosh$AC_PCA)
tirosh$AC_PCA_MICE = gsub('\\s+', '', tirosh$AC_PCA_MICE)
tirosh$OC_PCA = gsub('\\s+', '', tirosh$OC_PCA)
tirosh$OC_PCA_MICE = gsub('\\s+', '', tirosh$OC_PCA_MICE)

## Patient 
patients=c("MGH107","MGH135","MGH64","MGH45","MGH142_P1","MGH142_P2","MGH208","MGH208_P2","MGH201_P1","MGH201_P2")

for (m in 1:length(patients)) {
  ## Patient 
  patient = patients[m]
  print(patient)
  ## Read in TPM count matrices
  tpm = read.table( paste0( "/Gene_counts/",patient, ".tpm.txt")  )
  
  qc = read.table("IDH_filtered_cells.txt")
  rownames(qc)<-qc$V1
  
  tpm = tpm[ , grepl( paste( rownames(qc) , collapse = "|") ,  colnames(tpm)  )  ]
  
### Remove main genes from tpm then re-add after filtration
  mg = c()
  for (i in names(tirosh)) { mg = c(mg, as.character(tirosh[, i])[ complete.cases(as.character(tirosh[, i])) ]  ) }
  main_genes = tpm[ which(rownames(tpm) %in% mg) , ]
  tpm = tpm[ !(rownames(tpm) %in% rownames(main_genes)) , ]

### Filter out lowly expressed genes
tpm = tpm[ which( log2(rowMeans(tpm)+1) > 4 ), ]
## re-add the main genes
tpm <- rbind(tpm, main_genes)

## Aggregate express. for binning : E_a = log2( (tpm) + 1  ) 
Ea = log2(rowMeans(tpm)+1)

## Perform Eij = log2( (tpm/10) + 1  ) per Tirosh
tpm = log2( (tpm/10) + 1  )

## NOTE:  i x j matrix, where i = gene, j = cell, 'n' refers for cells in the equation
## Center expression levels, Er_ij = E_ij - avg( E_i..n )
tpm = sweep(tpm, 
            1 ,
            rowMeans(tpm),
            `-`)

## Define Control set of gene bins (30)
bins = bin(Ea, nbins = 30, method = "content")

## Module by Cell Data frame to hold cell scores
cell_scores = data.frame(matrix(NA, dim(tirosh)[2] , dim(tpm)[2] ))
colnames(cell_scores) = colnames(tpm)
rownames(cell_scores) = colnames(tirosh)

## Function to compute cell score given cell and module
SC = function( cell , geneset ){
  G_control = c()
  for (gene in geneset){
    genes_bin = names(Ea[ which(bins == bins[ which(names(Ea) == gene) ]) ])
    genes_bin_100 = if( !identical(genes_bin, character(0)) ) sample( genes_bin , size=100 ) else next
    genes_bin_expr = as.numeric( cell[ genes_bin_100 , ]  )
    G_control = c(G_control, genes_bin_expr)
  }
  cell_geneset_expr = mean( cell[geneset,] , na.rm = T)
  control_geneset_expr = mean(G_control) 
  return( cell_geneset_expr - control_geneset_expr )
}

for (i in colnames(tpm) ) {
  print(i)
  cell_name = i
  cell = dplyr::select(tpm, cell_name)
  for (j in colnames(tirosh) ){
    print(j)
    geneset = tirosh[ , j ][ complete.cases(tirosh[,j]) ]
    cell_scores[j,i] = SC(cell, geneset)
  }
}

master = cell_scores

## assign type
master = rbind(master, rep(NA, dim(master)[2]))
rownames(master)[dim(master)[1]] = "Type"

## Assign the type
for (cell in colnames(master) ){
  y = dplyr::select(master, cell)
  scores = sort( y[ -which( rownames(y) %in% c("Type") ) , ], decreasing = T)
  cell_type = rownames(y)[ which( y == scores[1] )  ]
  master["Type", cell ] = cell_type
}


master = data.frame(t(master))
master$cell = rownames(master)
master <- master[, c("cell","G1.S", "G2.M", "AC_PCA", "AC_PCA_MICE", "OC_PCA", "OC_PCA_MICE", "stemness" , "Type")]

master$new_stemness_PCA <- NA
for (i in 1:dim(master)[1]){
  master$new_stemness_PCA[i] <- as.numeric(as.character(master$stemness[i])) -  max(as.numeric(as.character(master$AC_PCA[i])), as.numeric(as.character(master$OC_PCA[i])))
}
## Assign new types including undifferentiated
master$new_Type_PCA = NA

for ( i in 1:dim(master)[1] ) {
  ### max of the 2 lineage scores and stemness  
  type_order <- c("AC_PCA",  "OC_PCA", "new_stemness_PCA")
  the_scores <- c(as.numeric(as.character(master$AC_PCA[i])), as.numeric(as.character(master$OC_PCA[i])),  as.numeric(as.character(master$new_stemness_PCA[i] )))
  max_score <- max(the_scores)
  second_highest <- max( the_scores[ !the_scores %in% max_score  ] )
  the_diff <- (max_score - second_highest)
  if( (max_score > 0.5) & (the_diff > 0.5 )   ) {
    master$new_Type_PCA[i] <- type_order[ which(the_scores == max_score)  ]
  } else {
    master$new_Type_PCA[i] <- "undifferentiated"
  }
}

### add new stemness score (MICE PCA)
master$new_stemness_PCA_MICE <- NA
for (i in 1:dim(master)[1]){
  master$new_stemness_PCA_MICE[i] <- as.numeric(as.character(master$stemness[i])) - max(as.numeric(as.character(master$AC_PCA_MICE[i])), as.numeric(as.character(master$OC_PCA_MICE[i])))
}
## Assign new types including undifferentiated
master$new_Type_PCA_MICE = NA

for ( i in 1:dim(master)[1] ) {
  ### max of the 2 lineage scores and stemness  
  type_order <- c("AC_PCA_MICE",  "OC_PCA_MICE", "new_stemness_PCA_MICE")
  the_scores <- c(as.numeric(as.character(master$AC_PCA_MICE[i])), as.numeric(as.character(master$OC_PCA_MICE[i])),  as.numeric(as.character(master$new_stemness_PCA_MICE[i] )))
  max_score <- max(the_scores)
  second_highest <- max( the_scores[ !the_scores %in% max_score  ] )
  the_diff <- (max_score - second_highest)
  if( (max_score > 0.5) & (the_diff > 0.5 )   ) {
    master$new_Type_PCA_MICE[i] <- type_order[ which(the_scores == max_score)  ]
  } else {
    master$new_Type_PCA_MICE[i] <- "undifferentiated"
  }
}

write.table(master, sprintf("%s_cellscores.txt", patients[m]), row.names = F, col.names = T, sep = "\t", quote = F  )

}
