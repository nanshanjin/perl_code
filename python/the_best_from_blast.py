import os,sys,re

def addtodict2(thedict, key_a, key_b, val): 
    if key_a in thedict:
        thedict[key_a].update({key_b: val})
    else:
        thedict.update({key_a:{key_b: val}})

mydata=file("At.vs.GbA.1e-5.m8.blast")
dicts={}
dicts2={}
for line in mydata:
    line=line.strip().split("\t")
    dicts2[line[0]+line[11]]=line
    #dicts.setdefault(line[0],[]).append(line[11])##see https://www.cnblogs.com/ywl925/p/3810598.html
    dicts.setdefault(line[0],{})[line[1]]=line[11]
#print dicts
bestfile=open("best.txt","w")
for key1 in dicts:
    #print dicts2[key1+max(dicts[key1])]
    #bestfile.write("%s\t%s\n" % (key1,max(dicts[key1])))
    #bestfile.write("\t".join(dicts2[key1+max(dicts[key1])])+"\n")
    maxkey=max(dicts[key1], key=dicts[key1].get)
    print key1,dicts[key1][maxkey]
