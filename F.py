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

def countFound(fromList,foundIn):
    count=0
    for fromElem in fromList:
        found=False
        for elem in foundIn:
            if(fromElem.interval.overlaps(elem.interval)):
                found=True
                break
        if(found): count+=1
    return count

if(len(sys.argv)!=3):
    exit(sys.argv[0]+" <predictions.bed> <gold-standard.bed>")
(predictionsFile,goldFile)=sys.argv[1:]

predictions=BedReader.readAll(predictionsFile)
gold=BedReader.readAll(goldFile)
numPredictions=len(predictions)
numGold=len(gold)
TP=countFound(predictions,gold)
FP=numPredictions-TP
FN=numGold-countFound(gold,predictions)
sens=float(countFound(gold,predictions))/float(numGold)
spec=float(countFound(predictions,gold))/float(numPredictions)
F=2*sens*spec/(sens+spec)
print("Sn="+str(sens)+" Sp="+str(spec)+" F1="+str(F))
print("TP="+str(TP)+" FP="+str(FP)+" FN="+str(FN))

