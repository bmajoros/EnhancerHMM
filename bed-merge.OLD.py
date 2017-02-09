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

if(len(sys.argv)!=2):
    exit(sys.argv[0]+" <in.bed>")
infile=sys.argv[1]

records=BedReader.readAll(infile)
records.sort(key=lambda x: x.interval.begin)
L=len(records)
i=0
while(i+1<L):
    thisRec=records[i]
    nextRec=records[i+1]
    if(thisRec.interval.overlaps(nextRec.interval)):
        thisRec.score=max(thisRec.score,nextRec.score)
        thisRec.interval.end=nextRec.interval.end
        del records[i+1]
        L-=1
        continue
    i+=1

for rec in records:
    print(rec.chr+"\t"+str(rec.interval.begin)+"\t"+str(rec.interval.end)
          +"\t"+rec.name+"\t"+str(rec.score))

