#!/usr/bin/env python
#=========================================================================
# This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
# License (GPL) version 3, as described at www.opensource.org.
# Copyright (C)2016 William H. Majoros (martiandna@gmail.com).
#=========================================================================
from __future__ import (absolute_import, division, print_function, 
   unicode_literals, generators, nested_scopes, with_statement)
from builtins import (bytes, dict, int, list, object, range, str, ascii,
   chr, hex, input, next, oct, open, pow, round, super, filter, map, zip)
# The above imports should allow this program to run in both Python 2 and
# Python 3.  You might need to update your version of module "future".
from Rex import Rex
rex=Rex()

MAX_P=0.00261067332642239 # FDR 0.05

#PATH="/gpfs/fs0/data/reddylab/Chris/script_BAC_analysis/cv19/CV19_sig"
PATH="/gpfs/fs0/data/reddylab/Chris/script_BAC_analysis/cv19/window1_CV19_data"

def findElements(array):
    elements=[]
    L=len(array)
    if(L==0): return elements
    begin=array[0][0]
    for i in range(1,L):
        if(array[i][0]!=array[i-1][0]+1):
            elements.append([begin,array[i-1][0]+1])
            begin=array[i][0]
    return elements

IN=open(PATH,"rt")
header=IN.readline()
bySubstrate={}
for line in IN:
    fields=line.split()
    (id,chr,start,stop,p,cov,dcov,ecov,fc,logp)=fields
    p=float(p)
    if(p>MAX_P): continue
    if(rex.find("\"(\S+)\"",chr)): chr=rex[1]
    array=bySubstrate.get(chr,None)
    if(not array): array=bySubstrate[chr]=[]
    array.append([int(start),p])
IN.close()

keys=bySubstrate.keys()
for chr in keys:
    array=bySubstrate[chr]
    array.sort(key=lambda x: x[0])
    elements=findElements(array)
    for elem in elements:
        print(chr+"\t"+str(elem[0])+"\t"+str(elem[1]))

