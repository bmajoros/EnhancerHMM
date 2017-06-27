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
import TempFilename

# Process command line:
if(len(sys.argv)!=2): exit(sys.argv[0]+" <in.fastb>")
(program,fastbFile)=sys.argv

# Read FASTB file:
fastb=Fastb(fastbFile)
numTracks=fastb.numTracks()
tempFiles=[]
for i in range(numTracks):
    track=fastb.getIthTrack(i)
    if(track.isDiscrete()): continue
    data=track.getData()
    name=track.getID()
    L=track.getLength()
    tmp=TempFilename.generate(".gnuplot")
    tempFiles.append([tmp,name])
    with open(tmp,"wt") as fh:
        for pos in range(L):
            fh.write(str(pos)+"\t"+str(data[pos])+"\n")

# Invoke gnuplot:
cmd='gnuplot -e "plot '
n=len(tempFiles)
for i in range(n):
    pair=tempFiles[i]
    file,name=pair
    cmd+='\\"'+file+'\\" with lines title \\"'+name+'\\"'
    if(i+1<n): cmd+=','
cmd+='; pause -1"'
os.system(cmd)

# Clean up:
for file,name in tempFiles:
    os.remove(file)

