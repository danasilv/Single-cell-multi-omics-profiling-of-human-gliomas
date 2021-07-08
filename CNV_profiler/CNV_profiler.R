suppressMessages(library(argparse))

### Arguments to script
parser <- ArgumentParser()
parser$add_argument( "--path" , type="character" , help = "Path to Cov Files")
parser$add_argument( "--tumor" , type="character" , help = "Prefix of Tumor Covs")
parser$add_argument( "--normal" , type="character" , help = "Prefix of Normal Covs")
parser$add_argument( "--chrom_sizes" , type="character" , help = "File with Chromosome Sizes")
parser$add_argument( "--resolution" , type="character" , help = "One of 'chromosome' or 'sliding'. [DEFAULT:chromosome]")
parser$add_argument( "--window_size" , type="character" , help = "If resolution = 'sliding', size of window (in MB) [DEFAULT=100]")
parser$add_argument( "--sliding_by" , type="character" , help = "If resolution = 'sliding', sliding by (in MB) [DEFAULT=100]")
parser$add_argument( "--output" , type="character" , help = "Prefix of Results")
parser$add_argument( "--region" , type="character" , help = "Return results of a specified region , Format: 'chr1:1200:1400' ")
args <- parser$parse_args()

## Load Packages
suppressMessages(library(data.table))
suppressMessages(library(GenomicRanges))
suppressMessages(library(dplyr))
suppressMessages(library(stringr))
specify_decimal <- function(x, k) trimws(format(round(x, k), nsmall=k))
as.numeric.factor <- function(x) {as.numeric(levels(x))[x]}

## Set Defaults
if( is.null(args$output) ) { args$output = args$tumor }
if( is.null(args$resolution) ) { args$resolution = 'chromosome' }
if( is.null(args$window_size) ) { args$window_size = 100 } else { args$window_size = as.numeric(args$window_size)   }
if( is.null(args$sliding_by) ) { args$sliding_by = 100 } else { args$sliding_by = as.numeric(args$sliding_by)  }

### Function to be used by script
find_region_cpgs <- function( covs, region){
  ## Turn covs into Granges 
  granges <- lapply(covs, function(sites) GRanges( seqnames = sites$chr, ranges = IRanges(start = sites$start, end = sites$end) , mcol = sites[ , c("meth", "c_reads", "t_reads") ] ) )
  ## Intersect Regions
  obj <- do.call(c, lapply( granges , function(y) { m <- as_tibble( subsetByOverlaps( y , region ) ) ; m <- m[ ,c("seqnames", "start", "end", "mcol.meth", "mcol.c_reads", "mcol.t_reads")] ; colnames(m)<- c("chr", "start", "end", "meth", "c_reads","t_reads") ; return(  length(m$meth)) } ) ) 
  ## Return the obj
  return(obj)
}
find_region_meth_cpgs <- function( covs, region){
  ## Turn covs into Granges 
  granges <- lapply(covs, function(sites) GRanges( seqnames = sites$chr, ranges = IRanges(start = sites$start, end = sites$end) , mcol = sites[ , c("meth", "c_reads", "t_reads") ] ) )
  ## Intersect Regions
  obj2 <- do.call(c, lapply( granges , function(y) { m <- as_tibble( subsetByOverlaps( y , region ) ) ; m <- m[ ,c("seqnames", "start", "end", "mcol.meth", "mcol.c_reads", "mcol.t_reads")] ; colnames(m)<- c("chr", "start", "end", "meth", "c_reads","t_reads") ; meth = length( which( m$meth == 100) ); return( meth)  } ) ) 
  return(obj2)
}
find_region_unmeth_cpgs <- function( covs, region){
  ## Turn covs into Granges 
  granges <- lapply(covs, function(sites) GRanges( seqnames = sites$chr, ranges = IRanges(start = sites$start, end = sites$end) , mcol = sites[ , c("meth", "c_reads", "t_reads") ] ) )
  ## Intersect Regions
  obj3 <- do.call(c, lapply( granges , function(y) { m <- as_tibble( subsetByOverlaps( y , region ) ) ; m <- m[ ,c("seqnames", "start", "end", "mcol.meth", "mcol.c_reads", "mcol.t_reads")] ; colnames(m)<- c("chr", "start", "end", "meth", "c_reads","t_reads") ; un_meth = length( which( m$meth == 0) ); return(  un_meth )  } ) ) 
  ## Return the obj
  return(obj3)
}


### Read in Arguments
path <- args$path
## Read in cov files
cov_files <- list.files(path = path , pattern = "*cov$")
cov_list <-  lapply(cov_files, function(x) fread( paste0(path,x) ) )
## Name columns of cov files
col_names = c("chr", "start", "end", "meth", "c_reads", "t_reads", "ID") 
cov_list = lapply(cov_list, setNames, col_names )
names(cov_list) <-  cov_files
## Convert list of cov dataframes to tibbles
tibble_files <- lapply( cov_list  , function(x) dplyr::as_tibble(x) )
covs <- tibble_files
### Set the Tumor Covs and Normal Covs
groups = c(args$tumor, args$normal)

### Separate Tumor from Normal
tumor_covs <- covs[ grep( paste0("^",groups[1]) , names(covs)) ]
normal_covs <- covs[ grep( paste0("^",groups[2]) , names(covs)) ]
### Remove tumor from normal if there are any

if (args$resolution == "chromosome")  {
  ## Perform the segmentation of the genome
  info <- read.table( args$chrom_sizes , header = T)
  ## regions
  regions <- info
  colnames(regions) <- c("chr", "start", "end")
  regions$chr <- paste0("chr", regions$chr)
}

### SLIDING WINDOW    

if ( args$resolution == "sliding" )  {

  if( is.null( args$region ) ) { 
    ## Perform the segmentation of the genome
    info <- read.table( args$chrom_sizes , header = T)
    ## regions
    regions <- c()
    for (i in seq(info$chr) ){
      tmp <- data.frame( seq(0, info$end[i], by = args$sliding_by*1e6 ) , seq(0, info$end[i], by = args$sliding_by*1e6) + args$window_size*1e6 )
      tmp <- cbind( rep(  paste0("chr",i)  , dim(tmp)[1]), tmp  )
      colnames(tmp) <- c("chr", "start", "end")
      regions <- rbind(regions, tmp)    
      }
  }
  
  if( !is.null( args$region ) ) { 
    #print(args$region)
    regions = data.frame( unlist(strsplit( args$region, ":"))[1], unlist(strsplit( args$region, ":"))[2], unlist(strsplit( args$region, ":"))[3] )
    colnames(regions) = c("chr", "start", "end")
    regions$start = gsub( ",", "", regions$start)
    regions$end = gsub( ",", "", regions$end)
    regions$start = as.numeric(regions$start)
    regions$end = as.numeric(regions$end)
    #print(regions)
    tmp = data.frame(   seq( regions$start, regions$end , by = args$sliding_by*1e6 )  ,  seq( regions$start, regions$end , by = args$sliding_by*1e6 ) + args$window_size*1e6 )
    tmp <- cbind( rep(  regions$chr , nrow(tmp) ), tmp  )
    colnames(tmp) = c("chr", "start", "end")
    regions = tmp
  }

}

### Perform the CNV analysis
results = data.frame( row.names = names(covs) )
results2 = data.frame( row.names = names(covs) )
results3 = data.frame( row.names = names(covs) )

for ( i in 1:dim(regions)[1] ){
  window <- regions[i, ]
  print(window)
  ## The obj below has cpgs / region/ cell
  new_obj <- find_region_cpgs( covs , GRanges( seqnames = window$chr , IRanges( window$start, window$end) ) )
  new_obj2 <- find_region_meth_cpgs( covs , GRanges( seqnames = window$chr , IRanges( window$start, window$end) ) )
  new_obj3 <- find_region_unmeth_cpgs( covs , GRanges( seqnames = window$chr , IRanges( window$start, window$end) ) )
  window_cpgs = data.frame(new_obj)
  window_meth_cpgs = data.frame(new_obj2)
  window_unmeth_cpgs = data.frame(new_obj3)

  colnames(window_cpgs) <- paste0( window$chr , ".", i)
  colnames(window_meth_cpgs) <- paste0( window$chr , ".", i)
  colnames(window_unmeth_cpgs) <- paste0( window$chr , ".", i)
  print(window_cpgs)
  print(window_meth_cpgs)
  print(window_unmeth_cpgs)
  results = cbind( results, window_cpgs)
  results2 = cbind(results2, window_meth_cpgs)
  results3 = cbind(results3, window_unmeth_cpgs)
  cat( sprintf( "%s%% Complete" , specify_decimal( 100*(i/dim(regions)[1]), k = 1 ) ) , "\n" ) 
}

## Write out number of cpgs
write.table( results , sprintf("%s_cpgs_w%s_s%s.txt", args$output , args$window_size, args$sliding_by ), col.names = T, row.names = T, quote = F, sep = "\t" )
write.table( results2 , sprintf("%s_cpgs_w%s_s%s_methCpGs.txt", args$output , args$window_size, args$sliding_by ), col.names = T, row.names = T, quote = F, sep = "\t" )
write.table( results3 , sprintf("%s_cpgs_w%s_s%s_unmethCpGs.txt", args$output , args$window_size, args$sliding_by ), col.names = T, row.names = T, quote = F, sep = "\t" )


