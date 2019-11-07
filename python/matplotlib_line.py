import re
import numpy as np
import pylab as pl
#import matplotlib.pyplot as plt

ref1=open('snp.num.in.repeat.area.per.50000bp.txt','r')
ref2=open('snp.num.in.cds.area.per.50000bp.txt','r')
ref4=open('gc.content.in.genome.per.50000bp.rate.up.txt','r')

######## figure snp num in repeat area
x1_axis=[]
y1_axis=[]
for line in ref1:
	line=line.strip('\n')
	a,b=re.split(r'\t',line)
	x1_axis.append(int(a))
	if b== '':
		b=0;
	y1_axis.append(int(b))
ref1.close()
######## figure snp num in cds area
x2_axis=[]
y2_axis=[]
for line in ref2:
	line=line.strip('\n')
	a,b=re.split(r'\t',line)
	x2_axis.append(int(a))
	if b== '':
		b=0;
	y2_axis.append(int(b))
ref2.close()
######## figure gc content num in genome
x4_axis=[]
y4_axis=[]
for line in ref4:
    line=line.strip('\n')
    a,b=re.split(r'\t',line)
    x4_axis.append(int(a))
    if b== '':
        b=0
    y4_axis.append(b)
ref4.close()

x1 = x1_axis# Make an array of x values
y1 = y1_axis# Make an array of y values for each x value
x2 = x2_axis
y2 = y2_axis
x4 = x4_axis
y4 = y4_axis
pl.figure(figsize=(12, 12))
plot1=pl.plot(x1, y1,'m',linewidth=1,label='snp number in repeat area')# use pylab to plot x and y
plot2=pl.plot(x2, y2,'g',linewidth=1,label='snp number in gene area')
plot4=pl.plot(x4, y4,'y',linewidth=1,label='gc content(*100) in genome')

pl.title('The distribution of repeat sequences, SNP in gene regions, SNP in repeat regions, and GC content in PXO genome',fontsize=8,color='k')# give plot a title
pl.xlabel('Pxo genome length(bp)',fontsize= 15,color= 'k')# make axis labels
pl.ylabel('number',fontsize= 15,color= 'k')
legend=pl.legend(shadow=True,fontsize='medium',loc='best') ##make legend
pl.xlim(0, 5300000)# set axis limits
pl.ylim(0, 450)
#pl.show()# show the plot on the screen
pl.savefig("examples.pdf")
pl.savefig("examples.png",dpi=800)
