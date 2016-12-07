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
import sys
import random
from SummaryStats import SummaryStats
from Rex import Rex
rex=Rex()

FACTOR=6
MIN_DIST=5
MAX_DIST=600
INCREMENT=1

def updateHash(fields,hash):
    if(not rex.find("(\S+)_(\d+)$",fields[0])):
        raise Exception("Can't parse first column: "+fields[0])
    chr=rex[1]
    pos=int(rex[2])
    if(not hash.get(chr,None)): hash[chr]=[]
    array=hash[chr]
    fields=list(map(lambda x:int(x),fields[1:]))
    array.append(fields)
    if(pos!=len(array)): exit(str(pos)+"!="+str(len(array)))

def computeCor(dist,hash):
    X=[]
    Y=[]
    for chr in hash.keys():
        array=hash[chr]
        L=len(array)
        last=L-dist
        for i in range(last):
            dnaThis=array[i][0]
            if(dnaThis<1): continue
            nextIndex=i+dist
            dnaNext=array[nextIndex][0]
            if(dnaNext<1): continue
            rnaThis=array[i][FACTOR]/dnaThis
            rnaNext=array[nextIndex][FACTOR]/dnaNext
            X.append(rnaThis)
            Y.append(rnaNext)
    return SummaryStats.correlation(X,Y)

#=========================================================================
#                                  main()
#=========================================================================
if(len(sys.argv)!=2):
    exit(sys.argv+" <edgeR-input.txt>")
infile=sys.argv[1]

# Read table
hash={}
IN=open(infile,"rt")
IN.readline() # header
while(True):
    line=IN.readline()
    if(not line): break
    fields=line.split()
    if(len(fields)!=8): continue
    updateHash(fields,hash)
IN.close()

# Consider different correlation lengths
for i in range(MIN_DIST,MAX_DIST,INCREMENT):
    c=computeCor(i,hash)
    print(i,c,sep="\t")


