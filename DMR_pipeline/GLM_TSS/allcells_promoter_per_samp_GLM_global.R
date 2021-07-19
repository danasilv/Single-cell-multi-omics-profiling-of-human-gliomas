library(dplyr)
library(ggplot2)
library(OneR )
library(BisRNA)
as.numeric.factor <- function(x) {as.numeric(levels(x))[x]}

### arg 1 = ends with ..promoters_by_cell_meth.rds 
### arg 2 = ends with ..promoters_by_cell_cpgs.rds
### arg 3 = group 1 (cell type or something else)
### arg 4 = group 2 

args = commandArgs(TRUE)

### Use this to perform DMR Analysis
meth = readRDS( args[1] )
meth$GENE <- as.character(meth$GENE)
meth <- meth[ !is.na(meth$GENE) ,]
## A gene by cell matrix (inputs are # of meth)
rownames(meth) <- meth$GENE
meth <- meth[, -1]

#### Use this to impose min CpG requirement
cpgs = readRDS( args[2] )
cpgs$GENE <- as.character(cpgs$GENE)
cpgs <- cpgs[ !is.na(cpgs$GENE) ,]
## A gene by cell matrix (inputs are # of CpGs)
rownames(cpgs) <- cpgs$GENE
cpgs <- cpgs[, -1]

## cell type vector
ct = c( args[3] ,  args[4] )

## Read in cell type information
cell_info <- read.table("GBM_cellscores.txt", header=T, stringsAsFactors = F )
cell_info$meth_cell <- gsub( ".cov", ".binarize.promoters.1kb.cov" , cell_info$meth_cell )

## do two groups comparison
cell_info$Type2[ grep("^AC|^MES", cell_info$Type1) ] <- "Diff"
cell_info$Type2[ grep("^NPC|^OPC", cell_info$Type1) ] <- "Stem"

## Take cells corresponding to cell types
chosen = cell_info

## Define Global pathway methylation (avg of all pathways per cell)
tmp_meth = mapply(as.numeric.factor, meth)
global = colMeans(tmp_meth, na.rm = T)

## tmp will hold your results
tmp = data.frame()

    for (x in seq(rownames(meth)) ){
                  print(x)
                  x = rownames(meth)[x]
       
                  x1 = meth[x,];
                  x2 = cpgs[x,];
                  
                  y1 = as.list( x1[ which( !is.na(x1) ) ] ) ; ## Remove NA's
                  y2 = as.list( x2[ names(y1) ] ) ; ## Remove NA's (num cpgs)
                  
                  if( length(y1) == 0 ) next
                  
                  ## Make DF
                  df <- data.frame( do.call(c, lapply( y1, as.numeric.factor)) , do.call(c, lapply( y2, as.numeric.factor))  ) ; 
                  colnames(df) <- c("meth", "cpgs");
                  
                  ## Add cell
		  df$cell = names(y1)
                  
		  ## assign groups
                  df$group <- NA
                  df$group <- chosen$Type2[ match(df$cell, chosen$meth_cell) ]                

	              df$group[ grep( sprintf("^%s", ct[1]) , df$group ) ] <- ct[1] ; 
                  df$group[ grep( sprintf("^%s", ct[2]) , df$group ) ] <- ct[2] ;

		  df <- df[ df$group %in% ct , ]
                  df <- df[ !is.na(df$group) , ]

		  if ( nrow(df) == 0 ) next

                  ## assign batch
                  df$batch = NA

                  ## cpgs per promoter filtration
                  df <- filter(df, cpgs >= 5)  
 
                  ## skip if no cells remain after your filtration
                  if ( nrow(df)==0 ) next
                  ## skip if only one cell type is represented after your filtration
                  if( length(unique(df$group)) != 2 ) next
                  ## skip if all of the pathway methylation values are 0
                  if ( all(df$meth == 0)  )  {next}
                  ## skipp if all of the pathway methylation values are 100
                  if ( all(df$meth == 1.0) ) {next}
                  
                  ## Minimum requirement = at least 10 cells in each group 
                  t1 = length(df$group[ df$group== ct[1] ])
                  t2 = length(df$group[ df$group== ct[2] ])
                 
		  if ( ! ( ( t1 >= 10 ) & ( t2 >= 10 ) ) ) next

                  ### Add global pathway methylation to the data frame
                  df$global <- global[ match(df$cell, names(global)) ]
                  
                  samp1 = df
                  if ( all(samp1$meth==0) ) next
                  ## b1 = group 1, b2 = group 2
                  b1 = samp1$meth[ which(samp1$group==ct[1]) ]
                  b2 = samp1$meth[ which(samp1$group==ct[2]) ]
                  ## Initial means
                  a1 = mean(b1)
                  a2 = mean(b2)
                  ## DELTAS
                  real_delta <- (a1-a2)
                  
                  ## Perform GLM
                  fit <- lm(meth ~ group + global , data = samp1)
                  res = summary(fit)$coefficients[,4]
                  p = res[ grep("group", names(res)) ]          
                  
                  #}
                  
                 
		  test = c( p, real_delta )

                  ans = data.frame( t(test) )
                  colnames(ans)[1] = "PVAL"
                  colnames(ans)[2] = "DELTA"
                  rownames(ans) = x
                  tmp = rbind( tmp, ans  )

              }
                  
analysis = sprintf("%s_%s_GLM_MINCPGS_5_GLOBAL", ct[1], ct[2])

## write out
tmp <- data.frame( rownames(tmp),  tmp)
colnames(tmp)[1] = "GENE"
rownames(tmp) = NULL
tmp = tmp[order(tmp$PVAL),]

### Save results
write.table(tmp, sprintf("allcells_DATA_TSS_%s.txt", analysis) , col.names = T, row.names = F, quote = F, sep = "\t")
