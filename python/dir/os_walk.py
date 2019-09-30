import os
import os.path
import sys
import re

def safe_open(file_name,mode):
    try:
        if not file_name.endswith('.gz'):
            return open(file_name,mode)
        else:
            import gzip
            return gzip.open(file_name,mode)
    except IOError:
        print file_name + ' do not exist!'
def create_dir (dir):
	if not os.path.exists(dir):
		assert not os.system('mkdir %s' % dir)

indir="/home/Account/yangns/project/metagenomics/JZ201907311146_Bacteria_Fungi/02-kraken2/"

header=["sampleID","Number_of_Bacterial_reads","Number_of_Fungal_reads","Fungi:Bacteria(%)"]
print "\t".join(header)
for parent,dirnames,filenames in os.walk(indir):
    for filename in filenames:
        if ".report.xls" in filename:
            files=filename.split(".")[0]
            filename=os.path.join(indir,filename)
            ba=0
            fu=0
            for line in safe_open(filename,"r"):
                if re.findall("Bacteria",line):
                    line=line.strip().split()
                    ba+=int(line[1])
                elif re.findall("Fungi",line):
                    line=line.strip().split()
                    fu+=int(line[1])
            fb="%.10f" % (float(fu)/ba)
            newline=[files]+[ba]+[fu]+[fb]
            print "\t".join(str(i) for i in newline)
