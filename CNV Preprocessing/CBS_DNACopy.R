library(DNAcopy)

#this code follows the instructions from the DNAcopy manual
dataA = read.csv('Chr7_CNV_for_R_DNAcopy_EGFR_25mb_downsampled.csv')

CNA.object <- CNA(cbind(dataA$MGH105A),
                  dataA$Chromosome,dataA$Position,
                  data.type="logratio",sampleid='MGH105A')

smoothed.CNA.object <- smooth.CNA(CNA.object)
segment.CNA.object <- segment(CNA.object, verbose=1,min.width=2)

output <- segment.CNA.object$output
p <- segments.p(segment.CNA.object, ngrid=100, tol=1e-6, alpha=0.05, search.range=100, nperm=1000)
write.csv(p,'Chr7_CBS_from_DNAcopy_EGFR_25mb_pvalues_MGH105A.csv')

write.csv(output,'Chr7_CBS_from_DNAcopy_DNAcopy_EGFR_25mb_MGH105A.csv')

pdf('Chr7_CBS_DNAcopy_DNAcopy_EGFR_25mb_MGH105A.pdf',width=6,height=4,paper='special')
plot(segment.CNA.object, plot.type="w")
dev.off()
