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
import os
import ProgramName
from BedReader import BedReader
from Rex import Rex
rex=Rex()

TIMEPOINTS=("t05","t1","t2","t3","t4","t5","t6","t7","t8","t10","t12")
BASE="/home/bmajoros/GGR/delta"
infile=BASE+"/genomewide-hg19.bed"

bySubstrate={}
reader=BedReader(infile)
while(True):
    record=reader.nextRecord() # Bed3Record or Bed6Record
    if(not record): break
    elem=None
    timepoint=None
    if(rex.find("^(\S+).standardized.*timepoints.(t\d+)",record.name)):
        elem=rex[1];
        timepoint=rex[2]
    elif(rex.find("^(\S+)\.(t\d+)\.fastb",record.name)):
        elem=rex[1];
        timepoint=rex[2]
    else: exit("cannot parse "+record.name)
    if(bySubstrate.get(elem,None) is None): bySubstrate[elem]=[]
    bySubstrate[elem].append([timepoint,record.score])
reader.close()

substrates=bySubstrate.keys()
for substrate in substrates:
    array=bySubstrate[substrate]
    numRecs=len(array)
    #if(numRecs!=12):
    if(numRecs<1):
        print(substrate,numRecs,sep="\t",file=sys.stderr)
        continue
    hash={}
    for rec in array:
        (timepoint,score)=rec
        hash[timepoint]=score
    print(substrate+"\t",end="")
    for timepoint in TIMEPOINTS:
        if(hash.get(timepoint,None) is None):
            exit(substrate+" has no "+timepoint)
        score=hash[timepoint]
        score=float(score)
        score=int(score)
        print(score,end="")
        if(timepoint!="t12"): print("\t",end="")
    print()
