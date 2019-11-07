#!/usr/bin/env python
# -*- coding: utf-8 -*-
u"""
File Name   : corr.py .

Author      : biolxy
E-mail      : biolxy@aliyun.com
Created Time: 2019-08-25 10:02:55
version     : 1.0
Function    : The author is too lazy to write nothing
Usage       :
"""
import os
import sys
import argparse
import re
import time
import json
import pandas as pd
from scipy import stats
import numpy as np
SCRIPT_FOLDER = os.path.abspath(os.path.dirname(__file__))
utilpath = os.path.dirname(SCRIPT_FOLDER)
sys.path.insert(1, utilpath)
from biolxyUtil.base import color_term
from biolxyUtil.base import mkdirf, progressPrintStr, MagicDict


__VERSION__ = 'v1.0'


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


def merge_lnc_mRNA_exp(gene_exp, lncRNA_exp, outputfile):
    OUT = open(outputfile, 'w')
    with open(gene_exp, 'r') as ff:
        for line in ff:
            if line.startswith("gene_id"):
                continue
            OUT.write(line)
    with open(lncRNA_exp, 'r') as ff:
        for line in ff:
            if line.startswith("transcript_id"):
                continue
            OUT.write(line)
    OUT.close()


def get_dfcorr2(df1, df2):
    """
    按行计算df1(lnc) 和 df2(gene) 中数据的相关性, 其中 df1、2的列数不许相等，method="pearson"
    100W 6 min 22.66 s
    """
    df1 = df1.loc[~(df1 == 0).all(axis=1), :]
    df2 = df2.loc[~(df2 == 0).all(axis=1), :]
    df1index = df1.index.tolist()
    df2index = df2.index.tolist()
    # list1, list2 = df1.values.tolist(), df2.values.tolist()
    # dict1 初始化
    # print(df1)
    dict1 = MagicDict()
    dict1['lnc'] = []
    dict1['gene'] = []
    dict1['corr'] = []
    dict1['pvalue'] = []
    for m in df1index:
        for n in df2index:
            corr, pvalue = stats.pearsonr(df1.loc[m], df2.loc[n])
            dict1['lnc'].append(m)
            dict1['gene'].append(n)
            dict1['corr'].append(corr)
            dict1['pvalue'].append(pvalue)
    df = pd.DataFrame.from_dict(dict1)
    return df


def get_node_edge(indf, indir, network_r, plotType, plotName):
    # plotType  mRNAlncRNA\ lncRNA\ circRNA\ circRNAmiRNA
    file1 = os.path.join(indir, "network.edge.txt")
    file2 = os.path.join(indir, "network.node.txt")
    df2 = indf[["lnc", "gene"]]
    df2["type"] = 1
    df2["cor"] = indf[["correlation"]]
    df2.rename(columns={'lnc': 'from', 'gene': 'to'}, inplace=True)
    df2.to_csv(file1, encoding='utf-8', sep='\t', index=False)
    df3 = indf[["gene"]]
    df3["media.type"] = 2
    df3["type.label"] = "mRNA"
    df3.rename(columns={'gene': 'id'}, inplace=True)
    df4 = indf[["lnc"]]
    df4["media.type"] = 1
    df4["type.label"] = "lncRNA"
    df4.rename(columns={'lnc': 'id'}, inplace=True)
    df4 = df4.append(df3)
    df5 = df4.drop_duplicates("id")
    df5.to_csv(file2, encoding='utf-8', sep='\t', index=False, mode='w')
    os.system(
        "Rscript {network_r} {edge} {node} {plotType} {plotName} ".format(
            network_r=network_r,
            edge=file1,
            node=file2,
            plotType=plotType,
            plotName=plotName))


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

def main():
    stepName = "correxp"
    echo = progressPrintStr(stepName)
    echo.start()

    # lncRNA_exp = sys.argv[1]
    # gene_exp = sys.argv[2]
    # 计算量为 1000 * 1000 = 100W
    # method 1 1m32s
    lncRNA_exp2 = os.path.join(outputdir, 'lncRNA_exp.f.xls')
    NumGreZero(lncRNA_exp, lncRNA_exp2, inpercent=0.5)
    df1 = pd.DataFrame(pd.read_csv(lncRNA_exp2, encoding='utf-8', sep='\t'))
    df1 = df1.set_index('transcript_id')
    gene_exp2 = os.path.join(outputdir, 'gene_exp.f.xls')
    NumGreZero(gene_exp, gene_exp2, inpercent=0.5)
    df2 = pd.DataFrame(pd.read_csv(gene_exp2, encoding='utf-8', sep='\t'))
    df2 = df2.set_index('gene_id')
    df3 = get_dfcorr(df1, df2)
    df3.to_csv(outputfile, encoding='utf-8', sep='\t', index=False)
    # method 2  7.6s
    # tmp = "tmp.xls"
    # merge_lnc_mRNA_exp(gene_exp, lncRNA_exp, tmp)
    # df = pd.DataFrame(pd.read_csv(tmp, encoding='utf-8', sep='\t', header=None))
    # df.set_index([0], inplace=True)
    # df2 = df.T
    # df3 = df2.corr()  # 格式还需要转换
    # method 3  16 min 7.01s
    # os.system(
    #     "Rscript /home/Account/lixy/scripts/rna-code/code/coexp-ceRNA/lincRNAmRNA_coexpression.R {lnc} {gene}"
    #     .format(lnc=lncRNA_exp, gene=gene_exp))
    # method 5  6min 20s
    # df1 = pd.DataFrame(pd.read_csv(lncRNA_exp, encoding='utf-8', sep='\t'))
    # df1 = df1.set_index('transcript_id')
    # df2 = pd.DataFrame(pd.read_csv(gene_exp, encoding='utf-8', sep='\t'))
    # df2 = df2.set_index('gene_id')
    # df3 = get_dfcorr2(df1, df2)
    df4 = df3[(df3.correlation >= corr) & (df3.pvalue <= pvalue)]
    df5 = df4.sort_values(by=['correlation'], ascending=False)  # 降序
    df5.to_csv(outputfile2, encoding='utf-8', sep='\t', index=False)
    indir = os.path.dirname(outputfile2)
    pdf = os.path.join(outputdir, "top400.diff_mRNA_lncRNA_co-exp")
    get_node_edge(df5.head(400), indir, network_r, "mRNAlncRNA", pdf)
    pdf = os.path.join(outputdir, "top500.diff_mRNA_lncRNA_co-exp")
    get_node_edge(df5.head(500), indir, network_r, "mRNAlncRNA", pdf)
    # end
    echo.stop()


if __name__ == '__main__':
    try:
        SCRIPT_FOLDER = os.path.abspath(os.path.dirname(__file__))
        parser = argparse.ArgumentParser(
            prog="correxp".format(__VERSION__),
            description=color_term(
                "This is a pipeline of sino-RNAseq, developed by biolxy and Yehao\
            in sinomics company, with all rights reserved.", "green", False))
        parser.add_argument('-r',
                            '--runoptions',
                            type=str,
                            help="Selection subroutine to run",
                            default='run',
                            metavar='')
        parser.add_argument(
            '-o',
            '--outputdir',
            type=str,
            help="the output file, such as diff_lncRNA_mRNA_cor_filiter.txt",
            metavar='')
        parser.add_argument('-l',
                            '--lncRNA',
                            type=str,
                            help="the expr file of lncRNA",
                            metavar='')
        parser.add_argument('-m',
                            '--mRNA',
                            type=str,
                            help="the expr file of mRNA",
                            metavar='')
        parser.add_argument('-c',
                            '--corr',
                            type=float,
                            help="the filiter num of corr",
                            metavar='')
        parser.add_argument('-p',
                            '--pvalue',
                            type=float,
                            help="the filiter num of pvalue",
                            metavar='')
        # parser.add_argument('-n',
        #                     '--number',
        #                     type=int,
        #                     help="the num of the network edge",
        #                     metavar='')
        args = parser.parse_args()
        lncRNA_exp = os.path.abspath(args.lncRNA)
        gene_exp = os.path.abspath(args.mRNA)
        outputdir = os.path.abspath(args.outputdir)
        outputfile = os.path.join(outputdir, "diff_lncRNA_mRNA_cor.txt")
        outputfile2 = os.path.join(outputdir,
                                   "diff_lncRNA_mRNA_cor_filiter.txt")
        corr = args.corr
        pvalue = args.pvalue
        network_r = os.path.join(SCRIPT_FOLDER, "network2.r")
        # number = args.number
        main()
    except KeyboardInterrupt:
        pass
    except IOError as e:
        raise
