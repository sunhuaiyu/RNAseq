# 20160509
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
# save('SummarizedExperiments', se)            
countTable = data.frame(assay(se))
rownames(countTable) = names(ebg)

library(Homo.sapiens)
countTable$symbol = mapIds(org.Hs.eg.db, keys=rownames(countTable), 
                           keytype='ENTREZID', column='SYMBOL')

# keep the rows with at least 5 reads in all samples together
write.csv(countTable[rowSums(countTable[, 2:11]) >=5, ],
           file='RNAseq_YuShi20160511.csv')
