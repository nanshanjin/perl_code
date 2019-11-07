
										  
#!/usr/bin/python
import sys,os
import argparse
from multiprocessing.dummy import Pool as ThreadPool

parser=argparse.ArgumentParser(description='split the big fastq file to small,default threadnum is 2')
parser.add_argument('-fq1','--fastq1',help='the fastq file',required=True)
parser.add_argument('-fq2','--fastq2',help='the fastq file',required=True)
parser.add_argument('-qcstat','--QCstatfile',help='the qc stat file of the fastq',required=True)
parser.add_argument('-n','--number',help='the number of split piece,default=2',type=int,default=2)
parser.add_argument('-o','--outdir',help='the out dir',required=True)
parser.add_argument('-p','--prefix',help='the prefix of split fastq',required=True)
argv=parser.parse_args()

'''
python split_fq.py
-fq1 MRDL1_WES_DHE02094-20_HJGJJCCXX_L3_1.clean.fq
-fq2 MRDL1_WES_DHE02094-20_HJGJJCCXX_L3_2.clean.fq
-o MRDL1_WES
-p MRDL1_WES_DHE02094-20_HJGJJCCXX_L3
-num 2

MRDL1_WES_DHE02094-20_HJGJJCCXX_L3-1_1.fq
MRDL1_WES_DHE02094-20_HJGJJCCXX_L3-1_2.fq
'''
if not os.path.exists(argv.fastq1):
	fastq1=(argv.fastq1) + ".gz"
	fastq2=(argv.fastq2) + ".gz"
else:
	fastq1=argv.fastq1
	fastq2=argv.fastq2

outdir=argv.outdir
qcstat=argv.QCstatfile
number=int(argv.number)
prefix=argv.prefix

if number==1:
	print 'When split number is 1,the script will do nothing!'
	exit()

def mv2map(fastq,outdir,qc_prefix):
	mapfastq=os.path.join(outdir,qc_prefix)
	if fastq.endswith('gz'):
		txt='gunzip -c %s>%s' %(fastq,mapfastq)
		assert not os.system('gunzip -c %s>%s' %(fastq,mapfastq))
	else:
		assert not os.system('ln -fs %s %s' %(fastq,mapfastq))
	return mapfastq

def extract_total_line(stat):
	stat=open(stat)
	stat.next()
	return stat.next().strip().split('\t')[-1]

def mv2newname(splitprefix,prefix,number,lane,outdir):
	a_name=[chr(i) for i in range(97,123)]
	for each in range(number):
		assert not os.renames('%sa%s' %(splitprefix,a_name[each]),'%s/%s-%d_%d.clean.fq' %(outdir,prefix,each+1,lane))

def split((fastq,outdir,qcstat,number,prefix,lane)):
	qc_prefix=os.path.basename(fastq).rstrip('.gz')
	mapfastq=mv2map(fastq,outdir,qc_prefix)
	N=int(extract_total_line(qcstat))
	n=((N/number)+1)*4
	splitprefix=mapfastq.lstrip('.fq')
	assert not os.system('split -l %d %s %s' %(n,mapfastq,splitprefix))
	assert not os.system('rm -rf %s' %(mapfastq))
	mv2newname(splitprefix,prefix,number,lane,outdir)
	


pool=ThreadPool(2)
pool.map(split,[(fastq1,outdir,qcstat,number,prefix,1),(fastq2,outdir,qcstat,number,prefix,2)])
pool.close()
pool.join()




