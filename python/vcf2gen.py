#!/home/Account/yangns/program/Python-2.7.14/bin/python 
# -*- coding" utf-8 -*-

import fileinput
import gzip
import string
import sys

if sys.version.startswith("3"):
    import io
    io_method = io.BytesIO
else:
    import cStringIO
    io_method = cStringIO.StringIO

import subprocess

sampleList=[]
limit1 = False


if len(sys.argv) > 1:
 aFile = sys.argv[1]
else:
 sys.exit()

#p = subprocess.Popen(["zcat", aFile], stdout = subprocess.PIPE)
#fh = io_method(p.communicate()[0])
#assert p.returncode == 0

fh = gzip.GzipFile(aFile, "r")

for line in fh:
#for line in fileinput.input(aFile, openhook=gzip.open):
  ss = string.split(line)

  if limit1:
     #['22', '16050435', '22:16050435', 'T', 'C', '.', 'PASS', 'MAF=0.00043;R2=0.01842', 'GT:DS:GP', 
     # '0/0:0.001:0.999,0.001,0.000', 
     chr = ss[0]
     srt = ss[1]
     nm  = ss[2]
     allel1 = ss[3]
     allel2 = ss[4]
     dosagesGT = ss[9:]

     print chr, nm, srt, allel1, allel2,
    
     for x in dosagesGT:
          y = x.split(':')
          print  string.join( y[2].split(',') ),
  
     print
     #break

  if "#CHROM" in ss[0]:
     sampleList = ss[9:]
#    for x in range(5, len(ss),3):
#      if ( float(ss[x+1]) >0 ):
#        print "SAMPLE ", (x-2)/3 , ":", "\t".join(ss[x:(x+3)])
     limit1 = True
     # write .sample list
     f = open('vcf.sample', 'w')
     for sample in sampleList:
        f.write("0 %s 0 0 0 -9\n" % sample)
     
     f.close()
