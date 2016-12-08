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
rex=Rex()

MIN_SIZE=10

def countElements(array,L,minSize):
    count=0
    begin=-1
    for i in range(L):
        if(i>0 and array[i-1]==0 and array[i]>0): begin=i
        elif(i==0 and array[i]>0): begin=i
        elif(i>0 and array[i-1]>0 and array[i]==0):
            length=i-begin
            if(length>=minSize):
                count+=1
                print(length)
            else:
                for j in range(begin,i): array[j]=0
    if(array[L-1]>0):
        length=L-begin
        if(length>=minSize):
            count+=1
            print(length)
        else:
            for j in range(begin,L): array[j]=0
    return count

def setHash(hash,chr,pos,value):
    if(hash.get(chr,None) is None): hash[chr]=[]
    array=hash[chr]
    if(array[pos] is None or value>array[pos]): array[pos]=value

def load(file,hash):
    with open(file,"rt") as fh:
        header=fh.readline()
        for line in fh:
            fields=line.split()
            if(len(fields)!=6): continue
            (label,fc,cpm,lr,p,fdr)=fields
            fc=float(fc); p=float(p); fdr=float(fdr)
            if(not rex.find("(\S+)_(\d+)$",label)):
                raise Exception("can't parse line: "+line)
            chr=rex[1]; pos=int(rex[2])-1
            if(fc>1 and fdr<=0.05):
                setHash(hash,chr,pos,fc)
            else: setHash(hash,chr,pos,0)

if(len(sys.argv)!=4):
    exit(sys.argv[0]+" <bac-lengths.txt> <in.eth> <in.dex>")
(lengthsFile,ethFile,dexFile)=sys.argv[1:]

hash={}
with open(lengthsFile,"rt") as IN:
    for line in IN:
        fields=line.split()
        if(len(fields)!=2): continue
        (chr,L)=fields
        hash[chr]=[0]*int(L)
load(ethFile,hash)
load(dexFile,hash)
totalElements=0
chroms=hash.keys()
for chr in chroms:
    array=hash[chr]
    L=len(array)
    totalElements+=countElements(array,L,MIN_SIZE)
    outfile=chr+".fastb"
    with open(outfile,"wt") as OUT:
        OUT.write("%logFC\n")
        for i in range(L):
            OUT.write(str(array[i])+"\n")
#print(totalElements,"elements found")
