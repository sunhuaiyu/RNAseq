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

# use cuffdiff instead
'''
cuffdiff -p 12 -o YS201606_cuffdiff/ -L YS01,YS02,YS03,YS04,YS05,YS07,YS08,YS10,YS11,YS13\
hg19_refgene.gtf \
YS01_sorted.bam YS02_sorted.bam YS03_sorted.bam YS04_sorted.bam \
YS05_sorted.bam YS07_sorted.bam YS08_sorted.bam YS10_sorted.bam \
YS11_sorted.bam YS13_sorted.bam >tmp.txt 2>tmp2.txt &
'''

# human gene names
hg = pd.read_csv('human_genes.txt', sep='\t', index_col='refseq')

fpkm = DataFrame({
     f[2:6]: pd.read_csv(f, sep='\t', index_col='gene_id')['FPKM'].sort_index()
     for f in glob('./YS*_cufflinks/genes.fpkm_tracking')})   

fpkm = pd.merge(fpkm, hg, left_index=True, right_index=True, how='left')
fpkm.index.name='RefSeq'
fpkm.to_csv('20160322YS_fpkm.csv')

'''
python fastq2FPKM_YS20160322.py > tmp.txt 2> tmp2.txt &
'''
# mapping statistics
f = glob('*Log.final.out')
for i in f:
    log = pd.read_csv(i, sep='\t', header=None)
    print i[:4], log.ix[8, 1], log.ix[7, 1], log.ix[4, 1]
'''
YS13 84.87% 24000419 28278560
YS11 84.79% 23470209 27681244
YS04 84.72% 18654020 22018002
YS05 85.19% 21519386 25259448
YS07 84.48% 21791774 25796596
YS01 85.38% 22649143 26528579
YS08 85.77% 17458565 20354415
YS02 85.20% 21635985 25394367
YS10 85.50% 22227623 25996128
YS03 84.91% 14955617 17613774
'''



###############################################
# differential expression YuShi RNAseq20160322
###############################################

import pandas as pd; from pandas import DataFrame, Series
import scipy.cluster.hierarchy as hier
from sklearn.preprocessing import scale

m = pd.read_csv('20160322YS_fpkm.csv', index_col='RefSeq')

#################################
# select for LIF responsive genes
# remove genes showing very low expression across all samples
select = m.ix[(m.ix[:, 0:4] >= 2).all(1), :] 
  #       & (m.ix[:, 0:4] >= 8).any(1), :] #this 2nd line is optional
#m = m.ix[(m.ix[:, 0:10] >= 0.01).sum(1) >= 6, :]; m.ix[:, 0:10] += 0.0001

select.ix[:, 'YS02_YS01'] = log2(select.YS02 / select.YS01)
select.ix[:, 'YS03_YS01'] = log2(select.YS03 / select.YS01)
select.ix[:, 'YS04_YS01'] = log2(select.YS04 / select.YS01)
LIFresponse = select.ix[
(abs(select.ix[:, ['YS02_YS01', 'YS03_YS01', 'YS04_YS01']]) >= 1).any(1), :]

# clustering and heatmap based on YS01 ~ YS04
data_matrix = scale(np.array(LIFresponse.ix[:, ['YS01', 'YS02', 'YS03', 'YS04']]), axis=1)

# hierarchical clustering
linkage_matrix = hier.linkage(data_matrix, method='ward')
heatmap_order = hier.leaves_list(linkage_matrix)
ordered_DataMatrix = data_matrix[heatmap_order,:]

# save a sorted FPKM table
LIFresponse = LIFresponse.ix[heatmap_order, :]
LIFresponse.to_csv('YS_selected_1a_LIFresponse_clustered.csv')

# draw heatmap by matplotlib
fig, ax = plt.subplots()
ax.imshow(ordered_DataMatrix, cmap=plt.cm.bwr, interpolation='nearest', aspect='auto')
ax.spines['left'].set_position(('outward', 10))
ax.spines['bottom'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.xaxis.set_ticks_position('none')
ax.yaxis.set_ticks_position('left')
plt.savefig('YS_selected_1a_LIFresponse_clustered.pdf')

#########################################
# select for IL6 responsive genes

select = m.ix[(m.ix[:, ['YS01', 'YS05', 'YS07']] >= 2).all(1), :]
#  & (m.ix[:, ['YS01', 'YS05', 'YS07']] >= 8).any(1), :] # this 2nd line is optional

select.ix[:, 'YS05_YS01'] = log2(select.YS05 / select.YS01)
select.ix[:, 'YS07_YS01'] = log2(select.YS07 / select.YS01)
IL6response = select.ix[(abs(select.ix[:, ['YS05_YS01', 'YS07_YS01']]) >= 1).any(1), :]

# hierarchical clustering based on YS01, YS05, YS07
data_matrix = scale(np.array(IL6response.ix[:, ['YS01', 'YS05', 'YS07']]), axis=1)
linkage_matrix = hier.linkage(data_matrix, method='ward')
heatmap_order = hier.leaves_list(linkage_matrix)
ordered_DataMatrix = data_matrix[heatmap_order,:]

# save a sorted FPKM table
IL6response = IL6response.ix[heatmap_order, :]
IL6response.to_csv('YS_selected_1c_IL6response_clustered.csv')

# draw heatmap by matplotlib
fig, ax = plt.subplots()
ax.imshow(ordered_DataMatrix, cmap=plt.cm.bwr, interpolation='nearest', aspect='auto')
ax.spines['left'].set_position(('outward', 10))
ax.spines['bottom'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.xaxis.set_ticks_position('none')
ax.yaxis.set_ticks_position('left')
plt.savefig('YS_selected_1c_IL6response_clustered.pdf')

##############################################

'''
# 1. download a program called "genePredToGtf" from here:
# http://hgdownload.cse.ucsc.edu/admin/exe/
# place the correct version of the executable somewhere in your path
 
# 2. Create the following file in your home directory:
echo 'db.host=genome-mysql.cse.ucsc.edu
db.user=genomep
db.password=password' > ~/.hg.conf
 
# the file's permissions must be user-only
chmod 0600 ~/.hg.conf
 
# 3. run "genePredToGtf" with any organism and any table that is in "genePred" format: 
## hg19/Ensemble Genes
genePredToGtf hg19 ensGene ensGene.gtf
  
# This will save "refGene.gtf" with all the required attributes 
# (gene_id, gene_name transcript_id, exon_number).
# but still not directly usable with DESeq because of multiple isoforms per gene.
"
'''

