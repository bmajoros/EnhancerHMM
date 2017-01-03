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

def readBySubstrate(filename):
    hash={}
    records=BedReader.readAll(filename)
    for rec in records:
        substrate=rec.chr
        if(hash.get(substrate,None) is None):
            hash[substrate]=[]
        hash[substrate].append(rec)
    return hash

if(len(sys.argv)!=3):
    exit(sys.argv[0]+" <sorted-predictions.bed> <sorted-gold-standard.bed>")
(predictionsFile,goldFile)=sys.argv[1:]

predictionsHash=readBySubstrate(predictionsFile)
goldHash=readBySubstrate(goldFile)

###
if(False):
    keys=goldHash.keys()
    for key in keys:
        if(predictionsHash.get(key,None) is None): continue
        predictions=predictionsHash[key]
        gold=goldHash[key]
        for prediction in predictions:
            category=find(prediction,gold)
            print(str(prediction.score)+"\t"+str(category))
    exit()
###

keys=goldHash.keys()
for key in keys:
    if(predictionsHash.get(key,None) is None): continue
    predictions=predictionsHash[key]
    gold=goldHash[key]
    numPred=len(predictions)
    numGold=len(gold)
    p=0
    g=0
    while(p<numPred and g<numGold):
        prediction=predictions[p]
        if(prediction.interval.overlaps(gold[g].interval)):
            print(str(prediction.score)+"\t1")
            p+=1
            #g+=1
            continue
        if(prediction.interval.end<=gold[g].interval.begin):
            print(str(prediction.score)+"\t0")
            p+=1
            continue
        if(gold[g].interval.end<=prediction.interval.begin):
            g+=1
            continue
        print(prediction.toString())
        print(gold[g].toString())
        raise Exception("internal error")
    while(p<numPred):
        prediction=predictions[p]
        print(str(prediction.score)+"\t0")
        p+=1
