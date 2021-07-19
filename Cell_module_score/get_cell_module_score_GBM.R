library(Seurat)
library(dplyr)
library(Homo.sapiens)
library(mygene)
library(OneR)

## Load in Neftel 2019 modules (GBM)
neftel = read.table("Neftel2019_metamodules.txt", sep = '\t', header = T, na.strings = "")
neftel$MES.1 = gsub('\\s+', '', neftel$MES.1)
neftel$MES.2 = gsub('\\s+', '', neftel$MES.2)
neftel$AC = gsub('\\s+', '', neftel$AC)
neftel$OPC = gsub('\\s+', '', neftel$OPC)
neftel$NPC.1 = gsub('\\s+', '', neftel$NPC.1)
neftel$NPC.2 = gsub('\\s+', '', neftel$NPC.2)
neftel$G1.S = gsub('\\s+', '', neftel$G1.S)
neftel$G2.M = gsub('\\s+', '', neftel$G2.M)

## ALL GBM PATIENTS
patients = c("MGH105A", "MGH105B", "MGH105C","MGH105D", "MGH121_P1", "MGH121_P2", "MGH121_P3", "MGH121_P4","MGH115", "MGH122", "MGH124","MGH211","MGH129")

for (m in 1:length(patients)) {
## Patient 
patient = patients[m]
print(patient)
## Read in TPM count matrices
tpm = read.table( paste0( "/Gene_counts/",patient, ".tpm.txt") )

qc = read.table("GBM_filtered_cells.txt")
rownames(qc)<-qc$V1

tpm = tpm[ , grepl( paste( rownames(qc) , collapse = "|") ,  colnames(tpm)  )  ]

### Remove main genes from tpm then re-add after filtration
  mg = c()
  for (i in names(neftel)) { mg = c(mg, as.character(neftel[, i])[ complete.cases(as.character(neftel[, i])) ]  ) }
  main_genes = tpm[ which(rownames(tpm) %in% mg) , ]
  tpm = tpm[ !(rownames(tpm) %in% rownames(main_genes)) , ]

### Filter out lowly expressed genes
tpm = tpm[ which( log2(rowMeans(tpm)+1) > 4 ), ]
## re-add the main genes
tpm <- rbind(tpm, main_genes)

## Aggregate express. for binning : E_a = log2( (tpm) + 1  ) 
Ea = log2(rowMeans(tpm)+1)

## Perform Eij = log2( (tpm/10) + 1  ) per Neftel
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
cell_scores = data.frame(matrix(NA, dim(neftel)[2] , dim(tpm)[2] ))
colnames(cell_scores) = colnames(tpm)
rownames(cell_scores) = colnames(neftel)

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
  for (j in colnames(neftel) ){
    print(j)
    geneset = neftel[ , j ]
    cell_scores[j,i] = SC(cell, geneset)
  }
}

master = cell_scores
master = rbind(master, rep(NA, dim(master)[2]))
master = rbind(master, rep(NA, dim(master)[2]))
rownames(master)[ dim(master)[1]-1  ] = "Type1"
rownames(master)[ dim(master)[1] ] = "Type2"

## Assign the type
for (cell in colnames(master) ){
  y = dplyr::select(master, cell)
  scores = sort(y[ -which( rownames(y) %in% c("Type1", "Type2", "G1.S", "G2.M") ) , ], decreasing = T)
  cell_type = rownames(y)[ which( y == scores[1] )  ]
  master["Type2", cell ] = cell_type
  
  y = dplyr::select(master, cell)
  scores = sort(y[ -which( rownames(y) %in% c("Type1", "Type2") ) , ], decreasing = T)
  cell_type = rownames(y)[ which( y == scores[1] )  ]
  master["Type1", cell ] = cell_type
  ## second type
}

master = data.frame(t(master))
master$cell = rownames(master)
write.table(master, sprintf("%s_cellscores.txt", patients[m]), row.names = F, col.names = T, sep = "\t", quote = F  )

}
