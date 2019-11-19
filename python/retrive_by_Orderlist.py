import re,sys,os
import pandas as pd

lists=[]
with open("/home/Account/cuilefang/prjt/1118/group_info_n.txt") as hand:
    for line in hand:
        if line.startswith("PATIENT_ID"):
            continue
        else:
            ids,ish=line.strip().split()
            lists.append(ids)
df = pd.read_table("/home/Account/cuilefang/prjt/1118/RSEM_normd.txt")
lists_filter=[]
lists_filter.append('Hugo_Symbol')
for line in lists:
    if line in list(df.columns):
        lists_filter.append(line)
#print list(df.columns)
#print lists_filter
#print df[lists_filter]
out=df[lists_filter]
out.to_csv("out.txt",sep="\t",index=0)
    
