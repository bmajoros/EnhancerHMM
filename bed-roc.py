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

def find(fromElem,foundIn):
    for elem in foundIn:
        if(fromElem.interval.overlaps(elem.interval)):
            return 1
    return 0

if(len(sys.argv)!=3):
    exit(sys.argv[0]+" <predictions.bed> <gold-standard.bed>")
(predictionsFile,goldFile)=sys.argv[1:]

predictions=BedReader.readAll(predictionsFile)
gold=BedReader.readAll(goldFile)
for prediction in predictions:
    category=find(prediction,gold)
    print(str(prediction.score)+"\t"+str(category))
