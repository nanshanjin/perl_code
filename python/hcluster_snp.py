##################################
########by nanshan.yang###########
##################################
import os 
import sys
import argparse
from argparse import RawTextHelpFormatter

parser = argparse.ArgumentParser(description="Hiraechical clustering samples with SNP.",formatter_class=RawTextHelpFormatter)
parser.add_argument('-i','--input',help="A list file contains germline snp, one sample per line, seperated by tab.\n   #sampleID\tsnp\n",required=True)
parser.add_argument('-o','--outdir',help="Output directory",required=True)
parser.add_argument('-r','--region',default='1',help="Chromosome or region(1:1-15000000) for hclust, default:1")
argv=parser.parse_args()


def safe_open(file_name,mode='r'):
	try:
		if not file_name.endswith('.gz'):
			return open(file_name,mode)
		else:
			import gzip
			return gzip.open(file_name,mode)
	except IOError:
		print file_name + ' do not exist!'

def retrieveVCF(vcfin,vcfout,chr,start,end):
	vcf_tmp = safe_open(vcfout,'w')
	if not os.path.exists(vcfin):
		vcfin = vcfin +'.gz'
		assert os.path.exists(vcfin)
	for line in safe_open(vcfin,'r'):
		if line.startswith('#'):
			vcf_tmp.write(line)
			continue
		array = line.strip().split('\t')
		pos = int(array[1])
		if argv.region and array[0] == chr and pos > start and pos < end:
			vcf_tmp.write(line)
	vcf_tmp.close()

if len(open(argv.input).read().strip().split('\n')) <= 2:
	sys.exit(1)

if not os.path.exists(argv.outdir):
	assert not os.system('mkdir '+argv.outdir)

region = [each.split('-') for each in argv.region.split(':')]
chr = region[0][0]
start = 0
end = 9999999999999
if len(region) > 1:
	start = int(region[1][0])
	end = int(region[1][1])
print chr
print start
print end
vcfdir = os.path.join(argv.outdir,'Vcfs')
if not os.path.exists(vcfdir):
	assert not os.system('mkdir '+vcfdir)


vcflist=[]
for line in open(argv.input):
	if line.startswith('#'):
		continue
	array = line.strip().split()
	vcfout = os.path.join(vcfdir,'%s.chr%s.snp.vcf' % (array[0],argv.region))
	retrieveVCF(array[1], vcfout, chr=chr,start=start,end=end)
	assert not os.system('bgzip -f %s' % vcfout)
	assert not os.system('tabix -p vcf %s.gz' % vcfout)
	vcflist.append(vcfout+'.gz')

mergedVCF = os.path.join(argv.outdir,'merged.chr%s.snp.vcf'%argv.region)
assert not os.system('vcf-merge %s \\\n>%s\n' % ('\\\n  '.join(vcflist), mergedVCF))

def gt2xls(gt):
	if '0/0' in gt:
		return 0
	elif '0/1' in gt:
		return 0.5
	elif '1/2' in gt:
		return 0.5
	elif '1/1' in gt:
		return 1
	elif gt == '.':
		return 0
	else:
		return 0
	

mergeXLS = os.path.join(argv.outdir,'merged.chr%s.snp.xls'%argv.region)
xls = safe_open(mergeXLS,'w')
samples = []
for line in safe_open(mergedVCF):
	if line.startswith('##'):
		continue
	array = line.strip().split('\t')
	if line.startswith('#CHROM'):
		samples = array[9:]
		xls.write('Point\t'+'\t'.join(samples)+'\n')
		continue
	xls.write('%s:%s\t'%(array[0],array[1])+'\t'.join([str(gt2xls(each)) for each in array[9:]])+'\n')
xls.close()

rscript = '''
setwd("%s")
dat<-read.table("merged.chr%s.snp.xls",head=T,sep="\\t",row.names=1)
dat<-data.matrix(dat)
d<-dist(t(dat))
h<-hclust(d)
png("Samples.cluster.png",width=min(ncol(dat)*50+200,3000),height=600,res=72*2,type='cairo-png')
plot(h,hang=-1,xlab="",ylab="",sub="")
dev.off()
'''%(argv.outdir, argv.region )
safe_open(os.path.join(argv.outdir,'hclust.R'),'w').write(rscript)
assert not os.system('Rscript %s' % (os.path.join(argv.outdir,'hclust.R')))
assert not os.system('rm -rf %s %s' % (vcfdir, mergedVCF))
