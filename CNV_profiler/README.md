README.md

Dependencies:
- R/3.6.1

Execution:
Rscript CNV_profiler.R --argument value 

Required arguments:
--path = "Path to .cov files from Bismark"
--tumor = "Prefix of Tumor .cov files"
--normal = "Prefix of Normal .cov files"
--chrom_sizes = "File with Chromosome Sizes"
--resolution = "One of 'chromosome' or 'sliding'. [DEFAULT:chromosome]"
--window_size = "If resolution = 'sliding', size of window (in MB) [DEFAULT=100]"
--sliding_by = "If resolution = 'sliding', sliding by (in MB) [DEFAULT=100]"
--output = "Prefix of Results"
--region = "Return results of a specified region , Format: 'chr1:1200:1400' "