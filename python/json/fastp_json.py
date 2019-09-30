import os
import sys
import re
import argparse
import json
# -*- coding: utf-8 -*-
from argparse import RawTextHelpFormatter

help='''
this script is to retrive NGS Q20\Q30 or others information from fastp json result file
'''
parser = argparse.ArgumentParser(description="this script is to retrive NGS Q20\Q30 or others information from fastp json result file\nContact:nanshan.yang@sinotechgenomics.com", formatter_class =RawTextHelpFormatter)
parser.add_argument('--pwd',help="Path of Project Directory for analysis(defalut: ./)",default=os.path.abspath('./'))
parser.add_argument('--jsons',help="The input json file by fastp\n",required=True)
parser.add_argument('--sample',help="The sample name \n",required=True)
argv = parser.parse_args()


def average(r1,r2,r3,r4):
    sums=r1+r2+r3+r4
    ave=sums/4
    return ave
def safe_open(file_name,mode='r'):
    try:
        if not file_name.endswith('.gz'):
            return open(file_name,mode)
        else:
            import gzip
            return gzip.open(file_name,mode)
    except IOError:
        print file_name + ' do not exist!'
jsonfile=file(argv.jsons)
raw_atgcn=safe_open(os.path.join(argv.pwd,argv.sample+".raw_atgcn"),"w")

#R1 raw
mydata=json.load(jsonfile)
r1_raw_a=mydata["read1_before_filtering"]["content_curves"]["A"]
r1_raw_t=mydata["read1_before_filtering"]["content_curves"]["T"]
r1_raw_c=mydata["read1_before_filtering"]["content_curves"]["C"]
r1_raw_g=mydata["read1_before_filtering"]["content_curves"]["G"]
r1_raw_n=mydata["read1_before_filtering"]["content_curves"]["N"]
i=0
##zip list
for r1_a,r1_t,r1_c,r1_g,r1_n in zip(r1_raw_a,r1_raw_t,r1_raw_c,r1_raw_g,r1_raw_n):
    #print r1_a,r1_t,r1_c,r1_g,r1_n
    raw_atgcn.write("%s\t%f\t%f\t%f\t%f\t%f\n" % (i,r1_a,r1_t,r1_c,r1_g,r1_n))
    i+=1
#R2 raw
r2_raw_a=mydata["read2_before_filtering"]["content_curves"]["A"]
r2_raw_t=mydata["read2_before_filtering"]["content_curves"]["T"]
r2_raw_c=mydata["read2_before_filtering"]["content_curves"]["C"]
r2_raw_g=mydata["read2_before_filtering"]["content_curves"]["G"]
r2_raw_n=mydata["read2_before_filtering"]["content_curves"]["N"]
j=0
for r2_a,r2_t,r2_c,r2_g,r2_n in zip(r2_raw_a,r2_raw_t,r2_raw_c,r2_raw_g,r2_raw_n):
    raw_atgcn.write("%s\t%f\t%f\t%f\t%f\t%f\n" % (j,r2_a,r2_t,r2_c,r2_g,r2_n))
    j+=1
raw_atgcn.close()
#R1 clean
clean_atgcn=safe_open(os.path.join(argv.pwd,argv.sample+".clean_atgcn"),"w")
r1_clean_a=mydata["read1_after_filtering"]["content_curves"]["A"]
r1_clean_t=mydata["read1_after_filtering"]["content_curves"]["T"]
r1_clean_c=mydata["read1_after_filtering"]["content_curves"]["C"]
r1_clean_g=mydata["read1_after_filtering"]["content_curves"]["G"]
r1_clean_n=mydata["read1_after_filtering"]["content_curves"]["N"]

#R2 clean
r2_clean_a=mydata["read2_after_filtering"]["content_curves"]["A"]
r2_clean_t=mydata["read2_after_filtering"]["content_curves"]["T"]
r2_clean_c=mydata["read2_after_filtering"]["content_curves"]["C"]
r2_clean_g=mydata["read2_after_filtering"]["content_curves"]["G"]
r2_clean_n=mydata["read2_after_filtering"]["content_curves"]["N"]
i2=0
for r1_a,r1_t,r1_c,r1_g,r1_n in zip(r1_clean_a,r1_clean_t,r1_clean_c,r1_clean_g,r1_clean_n):
    clean_atgcn.write("%s\t%f\t%f\t%f\t%f\t%f\n" % (i2,r1_a,r1_t,r1_c,r1_g,r1_n))
    i2+=1
j2=0
for r2_a,r2_t,r2_c,r2_g,r2_n in zip(r2_clean_a,r2_clean_t,r2_clean_c,r2_clean_g,r2_clean_n):
    clean_atgcn.write("%s\t%f\t%f\t%f\t%f\t%f\n" % (j2,r2_a,r2_t,r2_c,r2_g,r2_n))
    j2+=1
clean_atgcn.close()
#rawQ cleanQ
raw_Q=safe_open(os.path.join(argv.pwd,argv.sample+".raw_qual"),"w")
clean_Q=safe_open(os.path.join(argv.pwd,argv.sample+".clean_qual"),"w")
raw1q=mydata["read1_before_filtering"]["quality_curves"]["mean"]
raw2q=mydata["read2_before_filtering"]["quality_curves"]["mean"]
clean1q=mydata["read1_after_filtering"]["quality_curves"]["mean"]
clean2q=mydata["read2_after_filtering"]["quality_curves"]["mean"]
i3=0
for r1_q in raw1q:
    raw_Q.write("%s\t%f\n" % (i3,r1_q))
    i3+=1
i3=0
for r2_q in raw2q:
    raw_Q.write("%s\t%f\n" % (i3,r2_q))
    i3+=1
j3=0
for c1_q in clean1q:
    clean_Q.write("%s\t%f\n" % (j3,c1_q))
    j3+=1
j3=0
for c2_q in clean2q:
    clean_Q.write("%s\t%f\n" % (j3,c2_q))
    j3+=1
raw_Q.close()
clean_Q.close()
###stat
btotal_reads=mydata["summary"]["before_filtering"]["total_reads"]
btotal_bases=mydata["summary"]["before_filtering"]["total_bases"]
bq20_rate=mydata["summary"]["before_filtering"]["q20_rate"]
bq30_rate=mydata["summary"]["before_filtering"]["q30_rate"]
atotal_reads=mydata["summary"]["after_filtering"]["total_reads"]
atotal_bases=mydata["summary"]["after_filtering"]["total_bases"]
aq20_rate=mydata["summary"]["after_filtering"]["q20_rate"]
aq30_rate=mydata["summary"]["after_filtering"]["q30_rate"]
mystat=safe_open(os.path.join(argv.pwd,argv.sample+".stat"),"w")
mystat.write("\t".join(["#sampleID","rawreads","rawdata","rawq20","rawq30","cleanreads","cleandata","cleanq20","cleanq30"])+"\n")
mystat.write("\t".join([argv.sample,str(btotal_reads),str(btotal_bases),str(bq20_rate),str(bq30_rate),str(atotal_reads),str(atotal_bases),str(aq20_rate),str(aq30_rate)])+"\n")
mystat.close()
