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
from Fastb import Fastb
from FastbTrack import FastbTrack
from Rex import Rex
rex=Rex()

LOW_VALUE=3.5
HIGH_VALUE=4

#=========================================================================
# main()
#=========================================================================
if(len(sys.argv)!=5):
    exit(ProgramName.get()
         +" <in.bed> <label> <in-fastb-dir> <out-fastb-dir>\n")
(bedFile,label,inDir,outDir)=sys.argv[1:]

bedHash=BedReader.hashBySubstrate(bedFile)
files=os.listdir(inDir)
for file in files:
    id=file
    if(rex.find("(\S+).fastb",id)): id=rex[1]
    if(rex.find("(\S+)\.t\d+\.",id)): id=rex[1]
    elements=bedHash.get(id,[])
    fastb=Fastb(inDir+"/"+file)
    L=fastb.getLength()
    data=[LOW_VALUE]*L
    for elem in elements:
        for i in range(elem.interval.begin,elem.interval.end):
            data[i]=HIGH_VALUE
    track=FastbTrack("continuous",label,data,"")
    fastb.addTrack(track)
    fastb.save(outDir+"/"+file)
