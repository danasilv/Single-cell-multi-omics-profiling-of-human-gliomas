library(infercnv)

samplesGBM <- c("105A","105B","105C","105D","124","129","122","211","115","121")

for (sample in samplesGBM){
  wd = paste('/Users/lkluegel/Documents/GBM/scRNA/inferCNV_MGH',sample,sep='')
  input_raw = paste('MGH',sample,'.rsem.counts_inferCNV.txt',sep='')
  input_annotation = paste('MGH',sample,'_annotation_file_inferCNV.txt',sep='')
  input_order = paste('MGH',sample,'_gene_pos_inferCNV.txt',sep='')
  outd = paste('/Users/lkluegel/Documents/GBM/scRNA/inferCNV_MGH',sample,sep='')
  
  setwd(wd)
  
  infercnv_obj = CreateInfercnvObject(raw_counts_matrix=input_raw,
                                      annotations_file=input_annotation,
                                      delim="\t",
                                      gene_order_file=input_order,
                                      ref_group_names=c("normal"))
  
  infercnv_obj = infercnv::run(infercnv_obj,
                               cutoff=1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                               out_dir=outd,  # dir is auto-created for storing outputs
                               cluster_by_groups=T,   # cluster
                               denoise=T,
                               HMM=T
  )
}