# encoding: utf-8
import re,sys,os
import pandas as pd
import numpy as np
from scipy import stats
lnc= pd.read_table("Book_lncRNA.txt")
dicts={}
with open("list.tsv") as hand:
    for line in hand:
        if line.startswith("Samples"):
            continue
        else:
            s1,s2,g=line.strip().split()
            dicts.setdefault(g, []).append(s1)
with open("Combination.txt") as hand:
    for line in hand:
        if line.startswith("Condition1"):
            continue
        else:
            g1,g2,no=line.strip().split()
            newdata=lnc
            col_name=newdata.columns.tolist()
            #插入列
            col_name.insert(1,'Pvalue')
            newdata=newdata.reindex(columns=col_name)
            df1index = newdata[dicts[g1]].index.tolist()
            df2index = newdata[dicts[g2]].index.tolist()
            list1, list2 = newdata[dicts[g1]].values.tolist(), newdata[dicts[g2]].values.tolist()
            pv=[]
            for x,y in zip(list1,list2):
                a = np.array(x)
                b = np.array(y)
                #t检验的两组数据大小必须是一样的
                #https://blog.csdn.net/pipisorry/article/details/49515215
                #https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.ttest_rel.html#scipy.stats.ttest_rel
                s,p=stats.ttest_rel(a,b)
                #lnc['Pvalue']=p
                pv.append(p)
            newdata['Pvalue']=pv
            outfile=g1+"_vs_"+g2
            newdata.to_csv(outfile+".txt",sep="\t",index=0)
            newdata=lnc
