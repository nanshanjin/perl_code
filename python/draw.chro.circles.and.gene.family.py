### draw chromosome homology circles for eudicots with VV as reference, if VV is provided. If no vv chro is involved, it draws circles for other species
#coding=utf-8
import numpy as np
from matplotlib.patches import Ellipse, Circle
import matplotlib.pyplot as plt


import re
from math import *
from pylab import *
from matplotlib.patches import *
from bez import *

##### circle parameters
GAP_RATIO = 4 #gaps between chromosome circle, chr:gap = 4: 1
radius_a, radius_b = .33, .335    
sm_radius=(radius_b-radius_a)/2 #telomere capping
figurefile = "_".join(sys.argv[1:len(sys.argv)])+".genefam"

#### gene block parameters
blocklength = 0.01 # radian 0.01
block=0.01 #block scize
blockthick = 0.004 #0.006
shiftratio = -2.1 # define the distance between overlapping glocks

figure(1, (8, 8))  ### define the a square, or other rectangle of the figure, if to produce an oval here
root =axes([0, 0, 1, 1])

fpchrolen = file("Gm.lens")
fpgff = file("Gm.new.gff")
fpgenefamilyinf = file("Gm.ks.txt")
one=file("Gm.class.gene")

#print "Usage: python draw.chro.circles.py chrolist\n"

vvchro2color = {"VV1": "c", "VV2": "m", "VV3": "lightsalmon", "VV4": "aquamarine", "VV5": "orange", "VV6": "#FFFF33", "VV7": "#FFCC33","VV8": "#FF6666","VV9": "#FF3333","VV10": "#FF3399","VV11": "#FF99FF","VV12": "#CCFFCC","VV13": "#CCFF33","VV14": "#CCCCCC","VV15": "#CCCC99","VV16": "#CC99FF","VV17": "#CC9966","VV18": "#CC9900","VV19": "#CC66FF"}

#### specieslist and initial chrolist
chrolist = sys.argv[1:len(sys.argv)]
specieslist = []
iscompletegenome = {}
for i in range(len(chrolist)):
   string = chrolist[i]
#   print string[0:2]
   isnew = 1
   for sp in specieslist:
       if sp == string[0:2]:
       	isnew = 0
       break
   if isnew==1:
      specieslist.append(string[0:2])
      if string == string[0:2]:
         iscompletegenome[string[0:2]] = 1
      else:
           iscompletegenome[string[0:2]] = 0
#print "initial chro list", chrolist
#print "species list", specieslist
#print iscompletegenome

### input chromosome length
chro2len = {}
vvchrolist = []
otherchrolist = []
for row in fpchrolen:
    chro, length = row.split()
    chro = chro.upper()
    if(len(chro) > 5): continue
    sp = chro[:2].upper()
    if(iscompletegenome[sp] == 1):
                            chro2len[chro] = int(length)
                            if(re.search("VV", chro)):
                                               vvchrolist.append(chro)
                            else:
                                 otherchrolist.append(chro)
    else:
         if(chro in chrolist):
                 chro2len[chro] = int(length)
                 if(re.search("VV", chro)):
                                    vvchrolist.append(chro)
                 else:
                      otherchrolist.append(chro)
fpchrolen.close()
### full chro list
labels = vvchrolist + otherchrolist
#print "all chro", labels

## input gff
gene2pos={}
gene2chain = {}

for row in fpgff:
    ch, gene, start, end = row.split()
    start = int(start)
    end = int(end)
    gene2pos[gene] = int(start)
    gene2chain[gene] = int((end-start)/abs(end -start))
#    print gene, int(start), gene2chain
fpgff.close()

### input gene family gene1 gene2 Ka Ks
genepair2Ks = {}
genepair2Ka = {}
genes = []
for row in fpgenefamilyinf:
    gene1, gene2, Ka, Ks = row.split()
    genepair = ""
    if gene1 < gene2:
       genepair = gene1 + " " + gene2
    else:
       genepair = gene2 + " " + gene1
    if len(gene1.split("^")[0]) < 5 and gene1.split(".")[1] == "1" and len(gene2.split("^")[0]) < 5 and gene2.split(".")[1] == "1":
       genepair2Ka[genepair] = float(Ka)
       genepair2Ks[genepair] = float(Ks)
    if(gene1 not in genes):
	if len(gene1.split("^")[0]) < 5 and gene1.split(".")[1] == "1":
             genes.append(gene1)
    if(gene2 not in genes):
        if len(gene2.split("^")[0]) < 5 and gene2.split(".")[1] == "1":
             genes.append(gene2)
fpgenefamilyinf.close()

#print "genes ", genes, "genes"

def rad_to_coord(angle, radius):
    return radius*cos(angle), radius*sin(angle)

def to_deg(bp, total):
    # from basepair return as degree
    return bp*360./total

def to_radian(bp, total):
    # from basepair return as radian
#    print "to_radian", bp, total
    return radians(bp*360./total)

def plot_arc(start, stop, radius):
    # start, stop measured in radian
    #print start, stop
    t = arange(start, stop, pi/720.)
    x, y = radius*cos(t), radius*sin(t)
    plot(x, y, "k-", alpha=.5)

def plot_cap(angle, clockwise=True, radius=sm_radius):
    # angle measured in radian, clockwise is boolean
    if clockwise: t = arange(angle, angle+pi, pi/30.)
    else: t = arange(angle, angle-pi, -pi/30.)
    x, y = radius*cos(t), radius*sin(t)
    middle_r = (radius_a+radius_b)/2
    x, y = x + middle_r*cos(angle), y + middle_r*sin(angle)
    plot(x, y, "k-", alpha=.5)




fullchrolen = sum(chro2len.values())
chr_number = len(labels) # total number of chromosomes
GAP = fullchrolen/GAP_RATIO/chr_number # gap size in base pair
total_size = fullchrolen + chr_number * GAP # base pairs

start_list = [0]*chr_number
for i in xrange(1, chr_number):
    start_list[i] = start_list[i-1] + chro2len[labels[i-1]] + GAP
stop_list = [(start_list[i] + chro2len[labels[i]]) for i in xrange(chr_number)]


def transform_deg(ch, pos):
    return to_deg(pos + start_list[ch], total_size)

def transform_pt(ch, pos, r):
    # convert chromosome position to axis coords
#    print "transform", ch, pos, r
#    print "startlist", start_list[ch]
    rad = to_radian(pos + start_list[ch], total_size)
    return r*cos(rad), r*sin(rad)
def transform_pt2(rad, r):
    return r*cos(rad), r*sin(rad)

def plot_bez_inner(p1, p2, cl):
#    print "inner"
    a, b, c = p1
    ex1x, ex1y = transform_pt(a, b, c)
    a, b, c = p2
    ex2x, ex2y = transform_pt(a, b, c)
    # Bezier ratio, controls curve, lower ratio => closer to center
    ratio = .5
    x = [ex1x, ex1x*ratio, ex2x*ratio, ex2x]
    y = [ex1y, ex1y*ratio, ex2y*ratio, ex2y]
    step = .01
    t = arange(0, 1+step, step)
    xt = Bezier(x, t)
    yt = Bezier(y, t)
    plot(xt, yt, '-', color=cl, lw=.1, alpha=0.02)#alpha 

def plot_bez_outer(p1, p2, cl):
#    print "outer"
    a, b, c = p1
    ex1x, ex1y = transform_pt(a, b, c)
    a, b, c = p2
    ex2x, ex2y = transform_pt(a, b, c)
    # Bezier ratio, controls curve, lower ratio => closer to center
    ratio = 1.1
    x = [ex1x, ex1x*ratio, ex2x*ratio, ex2x]
    y = [ex1y, ex1y*ratio, ex2y*ratio, ex2y]
    step = .01
    t = arange(0, 1+step, step)
    xt = Bezier(x, t)
    yt = Bezier(y, t)
    plot(xt, yt, '-', color=cl, lw=.1)
def plot_arc_block(start, chain, radius,col):
    t = arange(start, start+block, pi/720.)
    x,y = radius * cos(t), radius*sin(t)
    x1, y1 = (radius-blockthick) * cos(t), (radius-blockthick) * sin(t)
    plot(x, y, "-", color=col, alpha=0.5, lw=1.5)
#    plot(x1, y1, "g-", alpha=0.5)

    x0, y0 = radius*cos(start), radius*sin(start)
    x1,y1  = (radius-blockthick)*cos(start), (radius-blockthick)*sin(start)
#    plot([x0, y0], [x1, y1], "g-", lw=0.2)

    x0, y0 = radius*cos(start+block), radius*sin(start+block)
    x1,y1  = (radius-blockthick)*cos(start+block), (radius-block)*sin(start+block)
#    plot([x0, y0], [x1, y1], "g-", lw=0.2)

def plot_bez_Ks(rad1, r1, rad2, r2, col):
#    print "bez Ks 1"
    ex1x, ex1y = transform_pt2(rad1, r1)
    ex2x, ex2y = transform_pt2(rad2, r2)
    ratio = 0.5#0.5
    x = [ex1x, ex1x*ratio, ex2x*ratio, ex2x]
    y = [ex1y, ex1y*ratio, ex2y*ratio, ex2y]
    step = .01
    t = arange(0, 1+step, step)
    xt = Bezier(x, t)
    yt = Bezier(y, t)
    plot(xt, yt, '-', color=col, lw=.1)

def plot_bez_Ks2(rad1, r1, rad2, r2, col):
#    print "bez Ks 2"
    ex1x, ex1y = transform_pt2(rad1, r1)
    ex2x, ex2y = transform_pt2(rad2, r2)
    ratio = -0.7#0.5
    sita = pi/2
    if ex1x != ex2x:
	sita = atan((ex2y-ex1y)/(ex2x-ex1x))
    d = sqrt((ex2x-ex1x)**2+(ex2y-ex1y)**2)
    L = d * ratio
    P1x = ex1x + L*sin(sita)
    P1y = ex1y - L*cos(sita)
    P2x = ex2x + L*sin(sita)
    P2y = ex2y - L*cos(sita)
    step = .01
    t = arange(0, 1+step, step)
    x=[ex1x, P1x, P2x, ex2x]
    y=[ex1y, P1y, P2y, ex2y]
    xt = Bezier(x,t)
    yt = Bezier(y,t)
    plot(xt, yt, '-', color = col, lw = 0.1)#0.1

def plot_bez_Ks3(rad1, r1, rad2, r2, col):
#    print "bez Ks 3"
    ex1x, ex1y = transform_pt2(rad1, r1)
    ex2x, ex2y = transform_pt2(rad2, r2)
    ratio = 0.5
    sita = pi/2
    if ex1x != ex2x:
	sita = atan((ex2y-ex1y)/(ex2x-ex1x))
    d = sqrt((ex2x-ex1x)**2+(ex2y-ex1y)**2)
    L = d * ratio
    P1x = ex1x + L*sin(sita)
    P1y = ex1y - L*cos(sita)
    P2x = ex2x + L*sin(sita)
    P2y = ex2y - L*cos(sita)
    step = .01
    t = arange(0, 1+step, step)
    x=[ex1x, P1x, P2x, ex2x]
    y=[ex1y, P1y, P2y, ex2y]
    xt = Bezier(x,t)
    yt = Bezier(y,t)
    plot(xt, yt, '-', color = col, lw = 0.1)
def plot_bez_Ks4(rad1, r1, rad2, r2, col):
#    print "bez Ks 3"
    ex1x, ex1y = transform_pt2(rad1, r1)
    ex2x, ex2y = transform_pt2(rad2, r2)
    ratio = -0.5
    sita = pi/2
    if ex1x != ex2x:
	sita = atan((ex2y-ex1y)/(ex2x-ex1x))
    d = sqrt((ex2x-ex1x)**2+(ex2y-ex1y)**2)
    L = d * ratio
    P1x = ex1x + L*sin(sita)
    P1y = ex1y - L*cos(sita)
    P2x = ex2x + L*sin(sita)
    P2y = ex2y - L*cos(sita)
    step = .01
    t = arange(0, 1+step, step)
    x=[ex1x, P1x, P2x, ex2x]
    y=[ex1y, P1y, P2y, ex2y]
    xt = Bezier(x,t)
    yt = Bezier(y,t)
    plot(xt, yt, '-', color = col, lw = 0.1)
def plot_bez_Ks5(rad1, r1, rad2, r2, col):
#    print "bez Ks 3"
    ex1x, ex1y = transform_pt2(rad1, r1)
    ex2x, ex2y = transform_pt2(rad2, r2)
    ratio = 0.4
    sita = pi/2
    if ex1x != ex2x:
	sita = atan((ex2y-ex1y)/(ex2x-ex1x))
    d = sqrt((ex2x-ex1x)**2+(ex2y-ex1y)**2)
    L = d * ratio
    P1x = ex1x + L*sin(sita)
    P1y = ex1y - L*cos(sita)
    P2x = ex2x + L*sin(sita)
    P2y = ex2y - L*cos(sita)
    step = .01
    t = arange(0, 1+step, step)
    x=[ex1x, P1x, P2x, ex2x]
    y=[ex1y, P1y, P2y, ex2y]
    xt = Bezier(x,t)
    yt = Bezier(y,t)
    plot(xt, yt, '-', color = col, lw = 0.1)
def plot_bez_Ks6(rad1, r1, rad2, r2, col):
#    print "bez Ks 3"
    ex1x, ex1y = transform_pt2(rad1, r1)
    ex2x, ex2y = transform_pt2(rad2, r2)
    ratio = -0.6
    sita = pi/2
    if ex1x != ex2x:
	sita = atan((ex2y-ex1y)/(ex2x-ex1x))
    d = sqrt((ex2x-ex1x)**2+(ex2y-ex1y)**2)
    L = d * ratio
    P1x = ex1x + L*sin(sita)
    P1y = ex1y - L*cos(sita)
    P2x = ex2x + L*sin(sita)
    P2y = ex2y - L*cos(sita)
    step = .01
    t = arange(0, 1+step, step)
    x=[ex1x, P1x, P2x, ex2x]
    y=[ex1y, P1y, P2y, ex2y]
    xt = Bezier(x,t)
    yt = Bezier(y,t)
    plot(xt, yt, '-', color = col, lw = 0.1)
def plot_sector(s):
    # block_id, chr_id, start_bp, stop_bp
    block_id, chr_id, start_bp, stop_bp = s
    theta0, theta1 = transform_deg(chr_id, start_bp), transform_deg(chr_id, stop_bp)
    dtheta = .1
    p1 = Wedge((0, 0), radius_a, theta0, theta1, dtheta=dtheta, fc="w", ec="w")
    p2 = Wedge((0, 0), radius_b, theta0, theta1, dtheta=dtheta, fc=colors[block_id], ec="w")
    root.add_patch(p2)
    root.add_patch(p1)


## sort gene according to lacation on circle
#
genessorted = []
geneno = 0
for i in xrange(len(genes)):
#    print i, genes[i]

    if geneno == 0:
       genessorted.append(genes[0])
       geneno = geneno + 1
    else:
         firstgene = genessorted[0]
         lastgene = genessorted[-1]
         chf = firstgene.split("^")[0]
         chl = lastgene.split("^")[0]
	 
         print firstgene, lastgene, chf, chl, gene2pos[firstgene]
         posf = gene2pos[firstgene] + start_list[labels.index(chf)]
         posl = gene2pos[lastgene] + start_list[labels.index(chl)]
         chi = genes[i].split("^")[0]
         posi = gene2pos[genes[i]] + start_list[labels.index(chi)]
#         print posf, posl, posi
         if posi <= posf:
            genessorted[0:0] = [genes[i]]
         elif posi >= posl:
              genessorted.append(genes[i])
         else:
           for j in xrange(len(genessorted)-1):
             chj = genessorted[j].split("^")[0]
             posj = gene2pos[genessorted[j]] + start_list[labels.index(chj)]
             chj1 = genessorted[j+1].split("^")[0]
             posj1 = gene2pos[genessorted[j+1]]+start_list[labels.index(chj1)]
             print posj, posj1, posi
             if posi > posj and posi < posj1:
                genessorted[j+1:j+1] = [genes[i]]
#    print "genesort ", genessorted

#    print "genesort ", genessorted


######
######
######the one gene
######
######

row = one.readline()#

oilgene={}
for row in one:
    onegene, fenlei = row.split()
    oilgene[onegene]=fenlei
    #rowno = rowno + 1
one.close()

######
######
######
######
######
# the chromosome layout
j = 0
for start, stop in zip(start_list, stop_list):
    start, stop = to_radian(start, total_size), to_radian(stop, total_size)
    # shaft
    plot_arc(start, stop, radius_a)
    plot_arc(start, stop, radius_b)

    # telemere capping
    plot_cap(start, clockwise=False)
    plot_cap(stop)
    
    # chromosome labels
    label_x, label_y = rad_to_coord((start+stop)/2, radius_b*0.9)#1.2
    #print label_x, label_y
    text(label_x, label_y, labels[j], horizontalalignment="center", verticalalignment="center", fontsize = 5, color = 'black')

    plot(label_x,label_y,'ro',markersize=15,lw=0.1,alpha=0.5, color="blue")#'ro'
    #plt.contour(label_x, label_y, label_x**1 + label_y**1, [1])
    #plt.axis('equal') 
    j+=1
########
########





### define shift level and draw gene blocks
shiftlevel = 0
gene2shift = {}
gene2location = {}
cho = genessorted[0].split("^")[0]
pos0 = gene2pos[genessorted[0]] + start_list[labels.index(cho)]
chain = gene2chain[genessorted[0]]
start = to_radian(pos0, total_size)
#plot_arc_block(start, chain, radius_a - shiftlevel * shiftratio * blockthick)
gene2location[genessorted[0]]= start
gene2shift[genessorted[0]] = shiftlevel
laststart = start

for i in xrange(len(genessorted)):
    chi = genessorted[i].split("^")[0]
    posi = gene2pos[genessorted[i]] + start_list[labels.index(chi)]
    chain = gene2chain[genessorted[i]]
    start = to_radian(posi, total_size)
    if(start-laststart < blocklength):
                       shiftlevel = shiftlevel + 1 
    #elif(0<start-laststart < 1):
    #     shiftlevel=shiftlevel + 1
    else:
        shiftlevel=0	 

    #shiftradius = radius_a - shiftlevel * shiftratio * blockthick

#    print "draw block: ", start, chain, shiftlevel, shiftratio, blockthick, shiftradius, radius_a

#    plot_arc_block(start, chain, radius_a - (shiftlevel+1) * shiftratio * blockthick, "blue")
    gene2location[genessorted[i]] = start
    gene2shift[genessorted[i]] = shiftlevel
    #print genessorted[i], chi, gene2pos[genessorted[i]], shiftlevel

    laststart = start


### Again: define shift level and draw gene blocks
shiftlevel = 0
gene2shift = {}
gene2location = {}
cho = genessorted[0].split("^")[0]
pos0 = gene2pos[genessorted[0]] + start_list[labels.index(cho)]
chain = gene2chain[genessorted[0]]
start = to_radian(pos0, total_size)
#plot_arc_block(start, chain, radius_a - shiftlevel * shiftratio * blockthick)
#gene2location[genessorted[0]]= start
#gene2shift[genessorted[0]] = shiftlevel
laststart = start

for i in xrange(len(genessorted)):
    chi = genessorted[i].split("^")[0]
    posi = gene2pos[genessorted[i]] + start_list[labels.index(chi)]
    chain = gene2chain[genessorted[i]]
    start = to_radian(posi, total_size)
    
    if(start-laststart<blocklength):
                       shiftlevel = shiftlevel + 1      
    else:
        shiftlevel=0	

    #shiftradius = radius_a - shiftlevel * shiftratio * blockthick
    #print "draw block: ", start, chain, shiftlevel, shiftratio, blockthick, radius_a
    if genessorted[i] in oilgene.keys():
        if(oilgene[genessorted[i]]=="2"):
            plot_arc_block(start, chain, radius_a - (shiftlevel+1) * shiftratio * blockthick,"darkorange")
        if(oilgene[genessorted[i]]=="1"):
            plot_arc_block(start, chain, radius_a - (shiftlevel+1) * shiftratio * blockthick,"darkblue")

 #   gene2location[genessorted[i]] = start
#    gene2shift[genessorted[i]] = shiftlevel
#    print genessorted[i], start

    laststart = start


#### draw tick showing scale of chromosomes
for i in xrange(chr_number):
   pos = 0
   while pos < chro2len[labels[i]]:
      xx1, yy1 = transform_pt(i, int(pos), radius_b)
      xx2, yy2 = transform_pt(i, int(pos), radius_b-0.003)
      plot([xx1, xx2], [yy1, yy2], "k-", lw = .2)
      pos = pos + 10000000

#



root.set_xlim(-.8, .8)#-.5, .5
root.set_ylim(-.8, .8)
root.set_axis_off()
savefig(figurefile + ".alpha.gene.pdf", dpi=1000)
savefig(figurefile + ".alpha.gene.png", dpi=1300)
