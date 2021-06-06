library(DNAcopy)

dataA = read.csv('Chr7_CNV_for_R_DNAcopy_EGFR_25mb_downsampled.csv')

CNA.object <- CNA(cbind(dataA$MGH105A),
                  dataA$Chromosome,dataA$Position,
                  data.type="logratio",sampleid='MGH105A')

smoothed.CNA.object <- smooth.CNA(CNA.object)
segment.CNA.object <- segment(CNA.object, verbose=1,min.width=2)

output <- segment.CNA.object$output
write.csv(output,'Chr7_CBS_from_DNAcopy_EGFR_25mb_MGH105A.csv')

pdf('Chr7_CBS_from_DNAcopy_EGFR_25mb_MGH105A.pdf',width=6,height=4,paper='special')
plot(segment.CNA.object, plot.type="w")
dev.off()

CNA.object <- CNA(cbind(dataA$MGH105B),
                  dataA$Chromosome,dataA$Position,
                  data.type="logratio",sampleid='MGH105B')

smoothed.CNA.object <- smooth.CNA(CNA.object)
segment.CNA.object <- segment(CNA.object, verbose=1,min.width=2)

output <- segment.CNA.object$output
write.csv(output,'Chr7_CBS_from_DNAcopy_EGFR_25mb_MGH105B.csv')

pdf('Chr7_CBS_from_DNAcopy_EGFR_25mb_MGH105B.pdf',width=6,height=4,paper='special')
plot(segment.CNA.object, plot.type="w")
dev.off()

CNA.object <- CNA(cbind(dataA$MGH105C),
                  dataA$Chromosome,dataA$Position,
                  data.type="logratio",sampleid='MGH105C')

smoothed.CNA.object <- smooth.CNA(CNA.object)
segment.CNA.object <- segment(CNA.object, verbose=1,min.width=2)

output <- segment.CNA.object$output
write.csv(output,'Chr7_CBS_from_DNAcopy_EGFR_25mb_MGH105C.csv')

pdf('Chr7_CBS_from_DNAcopy_EGFR_25mb_MGH105C.pdf',width=6,height=4,paper='special')
plot(segment.CNA.object, plot.type="w")
dev.off()

CNA.object <- CNA(cbind(dataA$MGH105D),
                  dataA$Chromosome,dataA$Position,
                  data.type="logratio",sampleid='MGH105D')

smoothed.CNA.object <- smooth.CNA(CNA.object)
segment.CNA.object <- segment(CNA.object, verbose=1,min.width=2)

output <- segment.CNA.object$output
write.csv(output,'Chr7_CBS_from_DNAcopy_EGFR_25mb_MGH105D.csv')

pdf('Chr7_CBS_from_DNAcopy_EGFR_25mb_MGH105D.pdf',width=6,height=4,paper='special')
plot(segment.CNA.object, plot.type="w")
dev.off()

CNA.object <- CNA(cbind(dataA$MGH115),
                  dataA$Chromosome,dataA$Position,
                  data.type="logratio",sampleid='MGH115')

smoothed.CNA.object <- smooth.CNA(CNA.object)
segment.CNA.object <- segment(CNA.object, verbose=1,min.width=2)

output <- segment.CNA.object$output
write.csv(output,'Chr7_CBS_from_DNAcopy_EGFR_25mb_MGH115.csv')

pdf('Chr7_CBS_from_DNAcopy_EGFR_25mb_MGH115.pdf',width=6,height=4,paper='special')
plot(segment.CNA.object, plot.type="w")
dev.off()

CNA.object <- CNA(cbind(dataA$MGH122),
                  dataA$Chromosome,dataA$Position,
                  data.type="logratio",sampleid='MGH122')

smoothed.CNA.object <- smooth.CNA(CNA.object)
segment.CNA.object <- segment(CNA.object, verbose=1,min.width=2)

output <- segment.CNA.object$output
write.csv(output,'Chr7_CBS_from_DNAcopy_EGFR_25mb_MGH122.csv')

pdf('Chr7_CBS_from_DNAcopy_EGFR_25mb_MGH122.pdf',width=6,height=4,paper='special')
plot(segment.CNA.object, plot.type="w")
dev.off()

CNA.object <- CNA(cbind(dataA$MGH124),
                  dataA$Chromosome,dataA$Position,
                  data.type="logratio",sampleid='MGH124')

smoothed.CNA.object <- smooth.CNA(CNA.object)
segment.CNA.object <- segment(CNA.object, verbose=1,min.width=2)

output <- segment.CNA.object$output
write.csv(output,'Chr7_CBS_from_DNAcopy_EGFR_25mb_MGH124.csv')

pdf('Chr7_CBS_from_DNAcopy_EGFR_25mb_MGH124.pdf',width=6,height=4,paper='special')
plot(segment.CNA.object, plot.type="w")
dev.off()

CNA.object <- CNA(cbind(dataA$MGH121),
                  dataA$Chromosome,dataA$Position,
                  data.type="logratio",sampleid='MGH121')

smoothed.CNA.object <- smooth.CNA(CNA.object)
segment.CNA.object <- segment(CNA.object, verbose=1,min.width=2)

output <- segment.CNA.object$output
write.csv(output,'Chr7_CBS_from_DNAcopy_EGFR_25mb_MGH121.csv')

pdf('Chr7_CBS_from_DNAcopy_EGFR_25mb_MGH121.pdf',width=6,height=4,paper='special')
plot(segment.CNA.object, plot.type="w")
dev.off()

CNA.object <- CNA(cbind(dataA$MGH129),
                  dataA$Chromosome,dataA$Position,
                  data.type="logratio",sampleid='MGH129')

smoothed.CNA.object <- smooth.CNA(CNA.object)
segment.CNA.object <- segment(CNA.object, verbose=1,min.width=2)

output <- segment.CNA.object$output
write.csv(output,'Chr7_CBS_from_DNAcopy_EGFR_25mb_MGH129.csv')

pdf('Chr7_CBS_from_DNAcopy_EGFR_25mb_MGH129.pdf',width=6,height=4,paper='special')
plot(segment.CNA.object, plot.type="w")
dev.off()

CNA.object <- CNA(cbind(dataA$MGH211),
                  dataA$Chromosome,dataA$Position,
                  data.type="logratio",sampleid='MGH211')

smoothed.CNA.object <- smooth.CNA(CNA.object)
segment.CNA.object <- segment(CNA.object, verbose=1,min.width=2)

output <- segment.CNA.object$output
write.csv(output,'Chr7_CBS_from_DNAcopy_EGFR_25mb_MGH211.csv')

pdf('Chr7_CBS_from_DNAcopy_EGFR_25mb_MGH211.pdf',width=6,height=4,paper='special')
plot(segment.CNA.object, plot.type="w")
dev.off()


