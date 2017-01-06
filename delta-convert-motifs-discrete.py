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
from Fastb import Fastb
from FastbTrack import FastbTrack

BASE="/home/bmajoros/GGR/delta"
MOTIFS_CONTINUOUS=BASE+"/motifs-continuous"
MOTIFS_DISCRETE=BASE+"/motifs-discrete"

def convert(infile,outfile):
    fastb=Fastb(infile)
    newTracks=[]
    numTracks=fastb.numTracks()
    for i in range(numTracks):
        oldTrack=fastb.getIthTrack(i)
        dataString=""
        for x in oldTrack.data:
            dataString+="1" if x>0 else "0"
        newTrack=FastbTrack("discrete",oldTrack.getID(),dataString,
                            oldTrack.deflineExtra);
        newTracks.append(newTrack)
    for i in range(numTracks):
        oldTrack=fastb.getIthTrack(0)
        fastb.dropTrack(oldTrack.getID())
    for newTrack in newTracks:
        fastb.addTrack(newTrack)
    fastb.save(outfile)

#=========================================================================
# main()
#=========================================================================
files=os.listdir(MOTIFS_CONTINUOUS)
for file in files:
    infile=MOTIFS_CONTINUOUS+"/"+file
    outfile=MOTIFS_DISCRETE+"/"+file
    convert(infile,outfile)



