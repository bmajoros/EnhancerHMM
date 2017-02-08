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
import sys
import os
import ProgramName
from Pipe import Pipe
from Rex import Rex
rex=Rex()

MUMMIE=os.environ["MUMMIE"]

def getParse(hmm,file):
    pipe=Pipe(MUMMIE+"/parse "+hmm+" "+file)
    pipe.readline() # discard header
    parse=[]
    while(True):
        line=pipe.readline()
        if(line is None): break
        state=int(line.rstrip())
        parse.append(state)
    return parse

def getPeaks(parse,peakStates):
    L=len(parse)
    begin=None
    peaks=[]
    for i in range(L):
        s=parse[i]
        if(s in peakStates):
            if(begin is None): begin=i
        else:
            if(begin is not None):
                peaks.append([begin,i])
                begin=None
    return peaks

def dumpPeaks(peaks,minPeakLen,substrate):
    for peak in peaks:
        (begin,end)=peak
        if(end-begin<minPeakLen): continue
        print(substrate,begin+1,end,sep="\t")

#=========================================================================
# main()
#=========================================================================
if(len(sys.argv)!=5):
    exit(ProgramName.get()
         +" <*.hmm> <peak-states> <min-peak-len> <fastb-dir>\n")
(hmmFile,peakStates,minPeakLen,fastbDir)=sys.argv[1:]
minPeakLen=int(minPeakLen)

# Parse list of peak states
fields=peakStates.rstrip().split(",")
peakStates=set()
for field in fields: peakStates.add(int(field))

# Proces FASTB directory
files=os.listdir(fastbDir)
for file in files:
    if(not rex.find("(\S+).fastb",file)): continue
    id=rex[1]
    parse=getParse(hmmFile,fastbDir+"/"+file)
    peaks=getPeaks(parse,peakStates)
    dumpPeaks(peaks,minPeakLen,id)

