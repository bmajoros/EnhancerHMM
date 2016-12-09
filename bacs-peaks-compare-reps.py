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
from Rex import Rex
from Interval import Interval
rex=Rex()

def countTrue(v):
    c=0
    for i in range(len(v)):
        if(v[i]): c+=1
    return c

def intersect(v1,v2):
    L=len(v1)
    V=[False]*L
    for i in range(L):
        V[i]=v1[i] and v2[i]
    return V

def union(v1,v2):
    L=len(v1)
    V=[False]*L
    for i in range(L):
        V[i]=v1[i] or v2[i]
    return V

def makeVector(peaks,chr,L):
    v=[False]*L
    for peak in peaks:
        if(peak.chr!=chr): continue
        for i in range(peak.begin,peak.end):
            v[i]=True
    return v

def loadPeaks(filename):
    array=[]
    with open(filename,"rt") as IN:
        for line in IN:
            fields=line.split()
            if(len(fields)<10): continue
            if(not rex.find("_peak_",fields[9])): continue
            (chr,start,end,length,abs_summit,pileup,pvalue,fold_enrichment,
             qvalue,name)=fields
            elem=Interval(int(start),int(end))
            elem.chr=chr
            elem.pvalue=pvalue
            elem.qvalue=qvalue
            elem.fold=fold_enrichment
            array.append(elem)
    array.sort(key=lambda i: i.begin)
    return array

#=========================================================================
# main()
#=========================================================================
if(len(sys.argv)!=4):
    exit(sys.argv[0]+" <bac-lengths.txt> <rep1.xls> <rep2.xls>")
(lengthsFile,rep1File,rep2File)=sys.argv[1:]

# Load bac lengths
lengths={}
with open(lengthsFile,"rt") as IN:
    for line in IN:
        fields=line.split()
        if(len(fields)!=2): continue
        (chr,length)=fields
        lengths[chr]=int(length)

# Load peaks
peaks1=loadPeaks(rep1File)
peaks2=loadPeaks(rep2File)
for chr in lengths.keys():
    v1=makeVector(peaks1,chr,lengths[chr])
    v2=makeVector(peaks2,chr,lengths[chr])
both=intersect(v1,v2)
either=union(v1,v2)
numer=countTrue(both)
denom=countTrue(either)
jaccard=float(numer)/float(denom)
print(str(jaccard)+" = "+str(numer)+"/"+str(denom))


