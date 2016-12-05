#!/bin/env python
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
import os
from Fastb import Fastb
from Pipe import Pipe
from Interval import Interval

BASE="/home/bmajoros/GGR/p300/"
HMM=BASE+"model/trained-pos5.hmm"
INDIR=BASE+"no-dna"
OUTDIR=BASE+"sliced"
FILE_LIST=BASE+"ls-raw.txt"
MUMMIE=os.environ["MUMMIE"]

#================================ FUNCTIONS ================================
def getStatePath(file):
    statePath=[]
    cmd=MUMMIE+"/parse "+HMM+" "+file
    pipe=Pipe(cmd)
    pipe.readline() # skip header
    while(True):
        line=pipe.readline()
        if(line is None): break
        state=int(line)
        statePath.append(state)
    return statePath

def getCoreEnhancer(path):
    L=len(path)
    pos=0
    while(pos<L):
        if(path[pos]!=1): break
        pos+=1
    begin=pos
    while(pos<L):
        if(path[pos]==5): break
        pos+=1
    end=pos
    return Interval(begin,end)

#================================ main() ================================
with open(FILE_LIST,"rt") as fh:
    for line in fh:
        file=line.rstrip()
        #if(file!="iter0_peak10082.standardized_across_all_timepoints.t10.fastb"): continue
        path=INDIR+"/"+file
        statePath=getStatePath(path)
        core=getCoreEnhancer(statePath)
        fastb=Fastb(path)
        fastb=fastb.slice(core.begin,core.end)
        fastb.save(OUTDIR+"/"+file)



