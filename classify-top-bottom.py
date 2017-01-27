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
import ProgramName
from Rex import Rex
rex=Rex()

if(len(sys.argv)!=3):
    exit(ProgramName.get()+" <genomewide.txt> <threshold>\n")
(infile,THRESHOLD)=sys.argv[1:]
THRESHOLD=float(THRESHOLD)

topCount={}
bottomCount={}
peaks=set()
IN=open(infile,"rt")
for line in IN:
    fields=line.rstrip().split()
    if(len(fields)!=7): continue
    (fastb,begin,end,id,LLR,parse,topBottom)=fields
    LLR=float(LLR)
    if(LLR<THRESHOLD): continue
    if(not rex.find("(\S+)\.t\d+\.fastb",fastb)): exit("can't parse "+fastb)
    peak=rex[1]
    hash=topCount if (topBottom=="top" or topBottom=="middle") else bottomCount
    hash[peak]=hash.get(peak,0)+1
    peaks.add(peak)
IN.close()
for peak in peaks:
    numTop=topCount.get(peak,0)
    numBottom=bottomCount.get(peak,0)
    #category="twopeak" if numTop>0 else "onepeak"
    category="twopeak" if numBottom>0 else "onepeak"
    print(peak,category,sep="\t")



