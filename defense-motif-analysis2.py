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
import ProgramName
from Interval import Interval
from Rex import Rex
rex=Rex()

MIN_PEAK_LEN=250
DELTA="/home/bmajoros/GGR/delta"
#ENHANCER_COORDS=DELTA+"/p300-hg38.bed"
PREDICTIONS=DELTA+"/genomewide-loopy.txt"
MOTIF_FASTBS=DELTA+"/motifs-continuous"

def loadPredictions(filename,wantTime):
    byChr={}
    IN=open(filename,"rt")
    for line in IN:
        fields=line.rstrip().split("\t")
        (peakID,begin,end,task,LLR,P,parse,topBottom)=fields
        if(not rex.find("(\S+)\.(t\d+)\.fastb",peakID)):
            raise Exception("can't parse"+peakID)
        peakID=rex[1]; thisTime=rex[2]
        if(thisTime!=wantTime): continue
        fields=parse.split("|")
        parse=[]
        for field in fields:
            if(not rex.find("(\d+):(\d+)-(\d+)",field)):
                raise Exception("can't parse "+field)
            state=int(rex[1]); begin=int(rex[2]); end=int(rex[3])
            if(state!=2 and state!=3 and state!=4): continue
            peak=Interval(begin,end)
            peak.type="peak" if state==3 else "hump"
            parse.append(peak)
        if(byChr.get(peakID,None) is None): byChr[peakID]=[]
        byChr[peakID].append(parse)
    IN.close()
    return byChr

#=========================================================================
# main()
#=========================================================================
if(len(sys.argv)!=2):
    exit(ProgramName.get()+" <timepoint>\n")
(timepoint,)=sys.argv[1:]

byChr=loadPredictions(PREDICTIONS,timepoint)



