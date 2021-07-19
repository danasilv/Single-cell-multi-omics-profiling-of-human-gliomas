library(data.table)
library(dplyr)

### List Cov Files
args = commandArgs(TRUE)
path = args[1]
files = list.files( path = path , pattern = "*bismark.cov$" )
files = paste0(path, files)

## Take small subet of files fo now
f = as.list(files)

## Given GENE, add reads counts per GENE 
add_reads <- function(x,m) { k = dplyr::filter( m , V11 == x ) ; val = paste0( sum(k$V5), ":", sum(k$V6) ) ; return( c(x, val) ) } 

covs <- lapply( f , function(x) { m <- fread(x) ; 
				  m$ID <- paste0( m$V1, m$V2, m$V3, m$V11) ; 
				  m <- m[ !duplicated(m$ID) , ]  ; 
				  ### Discard CpGs in between 0.9 and 0.1 meth
				  m <- filter(m, ( (V4 >= 90) | (V4 <= 10) ) )
				  new = m
				  ### more methylated reads = (1,0)                     
				  new$V5[  m$V5 > m$V6  ] <- 1
				  new$V6[  m$V5 > m$V6  ] <- 0				  
	              ### more unmethylated reads == (0,1)
				  new$V5[  m$V6 > m$V5  ] <- 0
                  new$V6[  m$V6 > m$V5  ] <- 1
				  write.table( new, gsub( ".cov", ".binarize.cov", x ) , col.names = F, row.names = F, quote = F, sep = "\t") 				  
				  ## out
				  #out <- sapply( unique( new$V11 ) , add_reads , new  ) ;
                                  #out <- data.frame( t(out) );
                                  #colnames(out) <- c("GENE", basename(x) ) ; 
				  return(NULL) } 
              )