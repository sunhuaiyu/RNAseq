# 2016-03-22 batch processing RNAseq fastq files
# need to prepare a table: sample name/code - fastq file name
import pandas as pd; from pandas import Series, DataFrame
import os
from glob import glob

sample = pd.read_csv('YuShi_RNAseq_files.txt', sep='\t', index_col='fastq_name')

link = 'http://igc2.snl.salk.edu/illumina/runs/160318-SE50/Hunter_Yu/'
os.system('wget -r -l1 -nd -nc -A.gz ' + link)
os.system('gunzip *.gz')

#STAR mapping
for f in glob('*.fastq'):
    os.system('STAR  --runThreadN 24  --genomeLoad LoadAndKeep' +
     '  --genomeDir /genomics/software/STAR/indexes/Genome_hg19' + 
     '  --readFilesIn ' + f +  
     '  --outFileNamePrefix ./' + sample.ix[f, 'code'] + '_STAR_' +    
     '  --outSAMstrandField intronMotif  --sjdbGTFfile ./hg19_refgene.gtf')

# sam -> bam -> sorted.bam file reformat
for f in glob('*.sam'):
    os.system('samtools view -b -S ' + f + ' > ' + f[:4] + '.bam' )

for f in glob('*.bam'):
    os.system('samtools sort -m 1000000000 ' + f + ' ' + f[:4] + '_sorted')    

# read counting with cufflinks
for f in glob('*_sorted.bam'):
    os.system('cufflinks -o ' + f[:4] + '_cufflinks/ -G ./hg19_refgene.gtf -p 24 ' + f)

# human gene names
hg = pd.read_csv('human_genes.txt', sep='\t', index_col='refseq')

fpkm = DataFrame({
     f[2:6]: pd.read_csv(f, sep='\t', index_col='gene_id')['FPKM'].sort_index()
     for f in glob('./YS*_cufflinks/genes.fpkm_tracking')})   

fpkm = pd.merge(fpkm, hg, left_index=True, right_index=True, how='left')
fpkm.index.name='RefSeq'
fpkm.to_csv('20160322YS_fpkm.csv')

'''
python RNAseq_fastq2FPKM_YS20160322.py > tmp.txt 2> tmp2.txt &
'''

