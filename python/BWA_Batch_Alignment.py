### BWA Batch Alignment ###

import os
from os import system

rawSeqFolder = 'raw_seq/' ## assume your raw data are stored in 'raw_seq' folder
bwaOutFolder = 'bwa_out/' ## assume the output would be in the 'bwa_out' folder

# change to your working directory which is one layer upper on rawSeqFolder and bwaOutFolder
os.chdir('~your working directory') 

seqFiles = os.listdir(rawSeqFolder)
seqPairs = {}

for fl in seqFiles:
    if not fl.endswith('.gz'):
        continue
    fileName = fl[:-11]
    if fileName not in seqPairs.keys():
        seq1 = rawSeqFolder+fileName+'_R1.fastq.gz'
        seq2 = rawSeqFolder+fileName+'_R2.fastq.gz'
        seqPairs[fileName] = seq1+' '+seq2
        
for pairName,seqData in seqPairs.items():
    bwaBam = bwaOutFolder+pairName+'.bam'
    sortedBam = bwaOutFolder+pairName+'.sorted.bam'
    os.system('bwa mem -t 8 -v 1 -M nipponbare7.fa '+seqData+'|samtools1.3 fixmate -O bam - '+bwaBam)
    os.system('samtools1.3 sort -@ 8 -T aligntmp -O bam -o '+sortedBam+' '+bwaBam) 

