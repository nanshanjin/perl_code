#!/usr/bin/python
import sys,os
import argparse
from multiprocessing.dummy import Pool as ThreadPool



def function((a,b,c,number,prefix,lane)):
    print a
	

pool=ThreadPool(2)
pool.map(function,[(a,b,c,number,prefix,1),(a,b,c,number,prefix,2)])
pool.close()
pool.join()

