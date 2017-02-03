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
from Fastb import Fastb
from FastbTrack import FastbTrack
from Rex import Rex
rex=Rex()

def getDnaTrack(peakID,dir):
    #fastb=Fastb(dir+"/"+peakID+".fastb")
    fastb=Fastb(dir+"/"+peakID+".standardized_across_all_timepoints.t00.fastb")
    return fastb.getTrackByName("dna")

#=========================================================================
# main()
#=========================================================================
if(len(sys.argv)!=4):
    exit(ProgramName.get()+" <in-fastb> <dna-dir> <out-fastb>\n")
(inDir,dnaDir,outDir)=sys.argv[1:]

files=os.listdir(inDir)
for file in files:
    if(not rex.find("(\S+)\.fastb",file)): exit("cant' parse filename: "+file)
    peakID=rex[1]
    dna=getDnaTrack(peakID,dnaDir)
    fastb=Fastb(inDir+"/"+file)
    fastb.addTrack(dna)
    fastb.save(outDir+"/"+file)



