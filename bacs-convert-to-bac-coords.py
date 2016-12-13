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
from BedReader import BedReader

if(len(sys.argv)!=4):
    exit(sys.argv[0]+" <in.bed> <bac-coords.txt> <out.bed>\n")
(inBed,bacCoordsFile,outFile)=sys.argv[1:]

# Load bac coords
bacs=[]
reader=BedReader(bacCoordsFile)
while(True):
    record=reader.nextRecord() # Bed3Record or Bed6Record
    if(not record): break
    bacs.append(record)
reader.close()

# Process input BED file
OUT=open(outFile,"wt")
reader=BedReader(inBed)
while(True):
    record=reader.nextRecord()
    if(not record): break
    for bac in bacs:
        if(bac.interval.contains(record.interval.begin)):
            record.interval.shift(-bac.interval.begin)
            record.chr=bac.name
            print(record.toString(),file=OUT)
reader.close()
OUT.close() 

       


