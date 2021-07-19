library(infercnv)

samplesGBM <- c("105A","105B","105C","105D","124","129","122","211","115","121")

for (sample in samplesGBM){
  wd = paste('/Users/lkluegel/Documents/GBM/scRNA/inferCNV_MGH',sample,sep='')
  
  setwd(wd)
  infercnv_obj = readRDS('run.final.infercnv_obj')
  
  write.csv(infercnv_obj@expr.data,file='preliminary_modified_expression.csv')
}
