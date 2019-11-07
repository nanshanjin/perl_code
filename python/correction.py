import os
import sys
import argparse
import re
import time
import pandas as pd
import numpy as np
def get_dfcorr(df1, df2):
    """
    按行计算df1(lnc) 和 df2(gene) 中数据的相关性, 其中 df1、2的列数不许相等，method="pearson"
    100W 1m32s
    """
    df1 = df1.loc[~(df1 == 0).all(axis=1), :]
    df2 = df2.loc[~(df2 == 0).all(axis=1), :]
    df1index = df1.index.tolist()
    df2index = df2.index.tolist()
    list1, list2 = df1.values.tolist(), df2.values.tolist()
    list_r = []
    for n, x in enumerate(list1):
        for m, y in enumerate(list2):
            a = np.array(x)
            b = np.array(y)
            r, p = stats.pearsonr(a, b)
            listtmp = [df1index[n], df2index[m], r, p]
            list_r.append(listtmp)
    df = pd.DataFrame(list_r, columns=["lnc", "gene", "correlation", "pvalue"])
    return df
def NumGreZero(infile, outfile, inpercent=0.5):
    OUT = open(outfile, 'w')
    inpercent = 0.9
    ff = open(infile, 'r')
    line = ff.readline()
    OUT.write("{}".format(line))
    for line in ff:
        line = line.rstrip()
        linelist = line.split("\t")
        numberlist = [float(x) for x in linelist[1:]]
        number = len([ x for x in numberlist if x > 0])
        if number * 1.0 / len(numberlist) > inpercent:
            OUT.write("{}\n".format(line))
    OUT.close()
lncRNA_exp2 = 'lncRNA_exp.f.xls'
NumGreZero('lncrnafile', lncRNA_exp2, inpercent=0.5)
df1 = pd.DataFrame(pd.read_csv(lncRNA_exp2, encoding='utf-8', sep='\t'))
df1 = df1.set_index('transcript_id')
gene_exp2 = 'gene_exp.f.xls'
NumGreZero('mrnafile', gene_exp2, inpercent=0.5)
df2 = pd.DataFrame(pd.read_csv(gene_exp2, encoding='utf-8', sep='\t'))
df2 = df2.set_index('gene_id')
df3 = get_dfcorr(df1, df2)
df3.to_csv(outputfile, encoding='utf-8', sep='\t', index=False)
