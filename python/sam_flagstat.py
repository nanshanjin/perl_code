#######################################################################
# Author           :	     nanshan yang  	
# Email            :         nanshangogo@163.com
# Filename         :		sam_flagstat.py
# Description      :		a_:array;d_:dict 	
#######################################################################
										  
#!/usr/bin/python
import os,sys
import HTSeq
import argparse
parser = argparse.ArgumentParser(description="samtools flagstat statistics")
parser.add_argument('--bam',help="the input bam file",type=str)
argv = vars(parser.parse_args())
bam = argv['bam']

stat = {'*':0,}

chrs = []
ta=np=dup=pe=se=nm=dchr=dchr_5=pm=0
for each in HTSeq.BAM_Reader(bam):
	if not each.iv:
		#stat['*'] += 1
		nm+=1
		continue
	chr = each.iv.chrom
	if chr not in stat:
		#stat[chr] = {'ta':0, 'np':0, 'dup':0, 'pe':0, 'se':0, 'na':0}
		chrs.append(chr)
	if each.not_primary_alignment:
		#stat[chr]['np'] += 1
		np+=1
		continue
	if each.proper_pair:
		pm+=1
	if each.aligned:
		#stat[chr]['ta'] += 1
		ta+=1
		if not each.mate_aligned:
			#stat[chr]['se'] += 1
			se+=1
		else:
			#stat[chr]['pe'] += 1
			pe+=1
			if chr != each.mate_start.chrom:
				dchr+=1
				if each.aQual >= 5:
					dchr_5+=1
	else:
		#stat[chr]['na'] += 1
		na+=1
	if each.pcr_or_optical_duplicate:
		#stat[chr]['dup'] += 1
		dup+=1
	
total=ta+nm
se=se+se
#print ta,np,dup,pe,se,nm,dchr,dchr_5
print 'Total:\t%d (100%%)' % total
print 'Duplicate:\t%d (%.2f%%)' % (dup, (dup*100.0)/total)
print 'Mapped:\t%d (%.2f%%)' % (ta, (ta*100.0)/total)
print 'Properly mapped:\t%d (%.2f%%)' % (pm, (pm*100.0)/total)
print 'PE mapped:\t%d (%.2f%%)' % (pe, (pe*100.0)/total)
print 'SE mapped:\t%d (%.2f%%)' % (se, (se*100.0)/total)
print 'With mate mapped to a different chr:\t%d (%.2f%%)' % (dchr, (dchr*100.0)/total)
print 'With mate mapped to a different chr ((mapQ>=5)):\t%d (%.2f%%)' % (dchr_5,(dchr_5*100.0)/total)

