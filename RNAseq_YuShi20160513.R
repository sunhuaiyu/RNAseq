# 20160513
# YuShi RNAseq sample analysis

sampleTable = read.csv('~/YuShi/YuShi_RNAseq_files.txt', sep='\t')
filenames = file.path('~/YuShi', paste0(sampleTable$code, '_sorted.bam'))

library("Rsamtools")
bamfiles = BamFileList(filenames)

library("GenomicFeatures")
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
ebg = exonsBy(TxDb.Hsapiens.UCSC.hg19.knownGene, by="gene")

library("GenomicAlignments")
library("BiocParallel")
register(SerialParam())

se = summarizeOverlaps(features=ebg, reads=bamfiles,
                        mode="Union",
                        singleEnd=FALSE,
                        ignore.strand=TRUE,
                        fragments=TRUE )
save(se, file='SummarizedExperiments20160513.RData')            
