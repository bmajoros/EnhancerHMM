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
from BedReader import BedReader
from Rex import Rex
from Fastb import Fastb
rex=Rex()

def scorePeak(peak,motifsDir):
    chr=peak.chr
    if(rex.find("(\S+).t\d+",chr)): chr=rex[1]
    motifFile=motifsDir+"/"+chr+".standardized_across_all_timepoints.t00.fastb"
    fastb=Fastb(motifFile)
    count=0
    for i in range(fastb.numTracks()):
        track=fastb.getIthTrack(i)
        intervals=track.getNonzeroRegions()
        for interval in intervals:
            if(interval.overlaps(peak.interval)): count+=1
    peakLen=peak.interval.length()
    score=float(count)/float(peakLen)
    return score

#=========================================================================
# main()
#=========================================================================
if(len(sys.argv)!=3):
    exit(ProgramName.get()+" <predictions.bed> <motifs-dir>\n")
(predictionsBed,motifsDir)=sys.argv[1:]

predictions=BedReader.readAll(predictionsBed)
for peak in predictions:

    ###
    peak.interval.begin=peak.interval.intCenter()-100
    peak.interval.end=peak.interval.intCenter()+100
    ###

    score=scorePeak(peak,motifsDir)
    print(score)



