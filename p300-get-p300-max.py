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
import os
import sys
import math
from BedReader import BedReader
from Pipe import Pipe
import SumLogProbs
from Interval import Interval
import ProgramName
import TempFilename
from Fastb import Fastb
from Rex import Rex
rex=Rex()

BASE="/home/bmajoros/GGR/p300"
fastbDir=BASE+"/no-dna"

def getMax(fastbFile):
    fastb=Fastb(fastbFile)
    track=fastb.getTrackByName("EP300")
    score=track.getMax()
    return score

#=========================================================================
# main()
#=========================================================================

if(len(sys.argv)!=2):
    exit(sys.argv[0]+" <partition.txt>")
infile=sys.argv[1]

IN=open(infile,"rt")
nextID=1
for line in IN:
    file=line.rstrip()
    if(not rex.find(".fastb$",file)): continue
    fullPath=fastbDir+"/"+file
    id="elem"+str(nextID)
    nextID+=1
    (p300,pos)=getMax(fullPath)
    begin=pos-200
    end=pos+200
    print(file,begin,end,id,p300,".",sep="\t")
IN.close()

