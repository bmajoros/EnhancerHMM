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

#=========================================================================
# main()
#=========================================================================
if(len(sys.argv)!=3):
    exit(ProgramName.get()+" <raw-dir> <out-dir>\n")
(rawDir,outDir)=sys.argv[1:]

files=os.listdir(rawDir)
for file in files:
    if(not rex.find("(iter0_peak\d+)\.standardized.*\.t00\.fastb",file)):
        continue
    peakID=rex[1]
    fastb=Fastb(rawDir+"/"+file)
    track=fastb.getTrackByName("dna")
    fastb=Fastb()
    fastb.addTrack(track)
    fastb.save(outDir+"/"+peakID+".fastb")



