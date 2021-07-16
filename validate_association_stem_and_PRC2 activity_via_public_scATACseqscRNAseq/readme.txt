The code is based on signac vignette from the satija lab. Here we re-analyzed data of GBM single-cell assay for transposase-accessible chromatin sequencing (scATACseq) from Wang et. al. GBM cells revealed clusters associated with the four core malignant cellular states described by scRNAseq in Neftel et. al. Gene expression activity inferred from scATAC-seq open chromatin (Methods) revealed a positive correlation between PRC2 targets accessibility and NPC/OPC-like cellular states in single-cells.

To run this code you will need 1. data folder containing scATACseq data (public from EGA, Wang et. al.) 2. gene signatures, including GBM signatures from Neftel et. al. and RC2 targets. These can be found in Neftel et. al. and Benporath et. al. and are also supplied in this github page under resources.

1. The folder should contain GBM signatures from the study: Neftel C, Laffy J, Filbin MG, Hara T, ... , Tirosh I, Suvà ML. An Integrative Model of Cellular States, Plasticity, and Genetics for Glioblastoma. Cell. 2019 Aug 8;178(4):835-849.e21. doi: 10.1016/j.cell.2019.06.024. Epub 2019 Jul 18. PMID: 31327527; PMCID: PMC6703186. A copy of these signatures is also in the github folder of this project: /validate_association_stem_and_PRC2 activity_via_public_scATACseqscRNAseq/resources/ 
resources.dir =  "..."

2. The folder should contain the publicly available data for tumor SF11956 (scATACseq: EGAF00002559211, scRNAseq: SF11956/RNA), from the study: Wang L, Babikir H, Müller S, Yagnik G et al. The Phenotypes of Proliferating Glioblastoma Cells Reside on a Single Axis of Variation. Cancer Discov 2019 Dec;9(12):1708-1719. PMID: 31554641
# To run it you will need to download or create the following files from EGA: 
# EGAF00002559211/aggr_normalized/outs/filtered_peak_bc_matrix.h5
# EGAF00002559211/aggr_normalized/outs/singlecell.csv
# EGAF00002559211/aggr_normalized/outs/fragments.tsv.gz
# EGAF00002559211/aggr_normalized/outs/fragments_filtered.tsv.gz 
# count tables can be found here: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE138794
# To project RNA signatures you will need the corresponfdind scRNAseq of tumor SF11956, available in EGA
# All data was created and published by Wang et. al. 
my.folder <- "..."
