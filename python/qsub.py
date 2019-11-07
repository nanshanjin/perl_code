#!/usr/bin/env python
# -*-coding:utf-8 -*-
"""
Contact:nanshan.yang@majorbio.com

Script:

Description:
this file is to qsub the shell script.
Only applicable for the shell script file.

Usage: this_program.py --sh input.sh

Options:
  --sh your input file 
  --version  Show version
  -h  help message
  
"""
from docopt import docopt
import sys,re,os,time
import argparse
from argparse import RawTextHelpFormatter
####
start=time.asctime()
####
parser = argparse.ArgumentParser(description="this file is the template python program.",formatter_class=RawTextHelpFormatter)
parser.add_argument('--sh',help="the shell script",required=True)
parser.add_argument('--resource',help="the resource like mem=3G",default='mem=3G')
parser.add_argument('--Queue',help="the Queue",default='dna')
parser.add_argument('--maxjobs',help="the number of maxjobs",default='20')
parser.add_argument('--CPU',help="set the cpu number to used for one job",default='4')
parser.add_argument('--nodes',help="set the procesor to submit for one jobs",default='1')
argv = vars(parser.parse_args())

def create_dir (dir):
	if not os.path.exists(dir):
		assert not os.system('mkdir %s' % dir)
def safe_open(file_name,mode='r'):
	try:
		if not file_name.endswith('.gz'):
			return open(file_name,mode)
		else:
			import gzip
			return gzip.open(file_name,mode)
	except IOError:
		print file_name + ' do not exist!'

#if not re.search('local',os.system('hostname')):
 #   sys.exit(1)

pid=os.getpid()#PID

if os.path.isfile(argv['sh']):
	shfile = argv['sh']
else:
	raise IOError('%s is not exsit' % argv['sh'])

#workdir=shfile+"."+pid+".qsub"

qsub=("qsub -l nodes=%s:ppn=%s -l %s -q %s %s" % (argv['nodes'],argv['CPU'],argv['resource'],argv['Queue'],shfile))
print qsub
logfile=file("logfile","w") % (pid)
logfile.write(qsub)
#safe_open(logfile,'w').write(qsub)
#os.system(qsub)

end=time.asctime()
print("\nOver!\n")
print ("start time: %s\n" % (start))
print("end time: %s\n" % (end))

