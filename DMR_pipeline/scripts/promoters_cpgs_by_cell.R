library(data.table)
library(dplyr)

### List Cov Files
args = commandArgs(TRUE)
path = args[1]
files = list.files( path = path , pattern = "*binarize.promoters.1kb.cov$" )
files = paste0(path, files)

f = as.list(files)

## Given GENE, add reads counts per GENE 
count_cpgs <- function(x,m) { k = dplyr::filter( m , V12 == x ) ; val = length(k$V5) ; return( c(x, val) ) } 

covs <- lapply( f , function(x) { m <- fread(x) ; 
				  m$ID <- paste0( m$V1, m$V2, m$V3, m$V12) ; 
				  m <- m[ !duplicated(m$ID) , ]  ; 
				  new = m
				  out <- sapply( unique( new$V12 ) , count_cpgs , new  ) ;
                                  out <- data.frame( t(out) );
                                  colnames(out) <- c("GENE", basename(x) ) ; 
				  return(out) } 
              )

fun = function(x,y) merge(x,y, by = c("GENE") , all = T) 

## Reduce cell
results <- Reduce(fun, covs)

saveRDS(results, sprintf("%s_binarize_promoters_cpgs.rds", args[2]) )