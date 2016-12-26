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
from Rex import Rex
rex=Rex()

BASE="/home/bmajoros/GGR/p300"
GENOMEWIDE=BASE+"/genomewide.txt"
P300_GR=BASE+"/p300-GR.bed"

p300GR=BedReader.readAll(P300_GR)
hasGR=set()
for rec in p300GR:
    hasGR.add(rec.name)

with open(GENOMEWIDE,"rt") as IN:
    for line in IN:
        fields=line.split()
        (fastb,begin,end,elemID,LLR,parse)=fields
        if(not rex.find("(iter0_peak\d+)",fastb)): exit("can't parse line")
        if(rex[1] in hasGR): print(line.rstrip())
        #if(rex[1] not in hasGR): print(line.rstrip())


