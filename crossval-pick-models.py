#!/usr/bin/env python
#=========================================================================
# This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
# License (GPL) version 3, as described at www.opensource.org.
# Copyright (C)2017 William H. Majoros (martiandna@gmail.com).
#=========================================================================
from __future__ import (absolute_import, division, print_function, 
   unicode_literals, generators, nested_scopes, with_statement)
from builtins import (bytes, dict, int, list, object, range, str, ascii,
   chr, hex, input, next, oct, open, pow, round, super, filter, map, zip)
# The above imports should allow this program to run in both Python 2 and
# Python 3.  You might need to update your version of module "future".
import os
from Rex import Rex
rex=Rex()

DELTA="/home/bmajoros/GGR/delta"
SLURMS=DELTA+"/slurms/crossval-slurms"
OUTPUTS=SLURMS+"/outputs"

def parseSlurm(fullPath):
    with open(fullPath,"rt") as IN:
        for line in IN:
            if(rex.find("crossval-pos-bin(\d+)-rep(\d+).hmm",line)):
                return [int(rex[1]),int(rex[2])]
    return None

def getLL(fullPath):
    best=None
    with open(fullPath,"rt") as IN:
        for line in IN:
            if(rex.find("\sLL=(\S+)",line)):
                LL=float(rex[1])
                if(best is None or LL>best): best=LL
    return best

#=========================================================================
# main()
#=========================================================================

likelihoods=[]
for i in range(5): likelihoods.append([])
slurms=os.listdir(SLURMS)
for slurm in slurms:
    if(not rex.find("command(\d+)",slurm)): continue
    commandIndex=rex[1]
    pair=parseSlurm(SLURMS+"/"+slurm)
    if(pair is None): continue
    (bin,rep)=pair
    outfile=OUTPUTS+"/"+commandIndex+".output"
    LL=getLL(outfile)
    likelihoods[bin-1].append([rep,LL])
for i in range(5):
    bin=likelihoods[i]
    bin.sort(key=lambda x: -x[1])
    for elem in bin: print(elem[0],elem[1])
    bestPair=bin[0]
    (rep,LL)=bestPair
    binNumber=i+1
    print(binNumber,rep,LL,sep="\t")
        
