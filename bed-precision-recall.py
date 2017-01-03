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

def readBySubstrate(filename):
    hash={}
    records=BedReader.readAll(filename)
    for rec in records:
        substrate=rec.chr
        if(hash.get(substrate,None) is None):
            hash[substrate]=[]
        hash[substrate].append(rec)
    #for key in hash.keys():
    #    hash[key].sort(key=lambda x: x.score)
    return hash

def getAllScores(hash):
    scores=set()
    keys=hash.keys()
    for key in keys:
        array=hash[key]
        for elem in array:
            scores.add(elem.score)
    array=list(scores)
    array.sort()
    return array

def getThresholds(scores):
    thresholds=[]
    N=len(scores)
    if(N<1): exit("no elements")
    thresholds.append(scores[0]-1.0)
    for i in range(N-1):
        thisScore=scores[i]
        nextScore=scores[i+1]
        threshold=(thisScore+nextScore)/2.0
        thresholds.append(threshold)
    thresholds.append(scores[N-1]+1.0)
    return thresholds

def filter(predictions,threshold):
    filtered=[]
    numPred=len(predictions)
    for i in range(numPred):
        prediction=predictions[i]
        if(prediction.score>=threshold):
            filtered.append(prediction)
    return filtered

#=========================================================================
# main()
#=========================================================================
if(len(sys.argv)!=3):
    exit(sys.argv[0]+" <sorted-predictions.bed> <sorted-gold-standard.bed>")
(predictionsFile,goldFile)=sys.argv[1:]

# Read files
predictionsHash=readBySubstrate(predictionsFile)
goldHash=readBySubstrate(goldFile)
chroms=list(goldHash.keys())
#chroms.sort()

# Get thresholds
scores=getAllScores(predictionsHash)
thresholds=getThresholds(scores)

# Iterate over thresholds
for threshold in thresholds:
    TP=0; FP=0; FN=0; totalPredictions=0; totalGold=0
    for chr in chroms:
        if(predictionsHash.get(chr,None) is None): continue
        predictions=filter(predictionsHash[chr],threshold)
        gold=goldHash[chr]
        numPredictions=len(predictions)
        totalPredictions+=numPredictions
        numGold=len(gold)
        totalGold+=numGold
        tp=countFound(predictions,gold)
        fp=numPredictions-tp
        fn=numGold-countFound(gold,predictions)
        TP+=tp; FP+=fp; FN+=fn
    recall=TP/float(totalGold)
    precision=TP/float(totalPredictions) if\
        totalPredictions>0 else 0.0
    print(recall,precision,sep="\t")
