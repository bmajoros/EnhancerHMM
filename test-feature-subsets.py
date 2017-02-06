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
from Pipe import Pipe

# Global constants
NUM_PARTITIONS=5
DELTA="/home/bmajoros/GGR/delta"
HMM_DIR=DELTA+"/hmm/crossval/best"
BEST_MODELS=(210,110,33,214,249)
TRAIN_POS=DELTA+"/fastb/crossval/training-sets/pos"
TRAIN_NEG=DELTA+"/fastb/crossval/training-sets/neg"
TEST_POS=DELTA+"/fastb/crossval/partitions/pos"
TEST_NEG=DELTA+"/fastb/crossval/partitions/neg"
TEMP_POS_HMM=DELTA+"/hmm/subset-pos.hmm"
TEMP_POS_HMM=DELTA+"/hmm/subset-neg.hmm"
TEMP_POS_TEST=DELTA+"/subset-pos"
TEMP_NEG_TEST=DELTA+"/subset-neg"
MUMMIE=os.environ["MUMMIE"]

# Feature sets
ALL_FEATURES=["DNase.t00", "H3K27ac.t00", "H3K4me1.t00", "H3K4me2.t00", 
              "EP300.t00", "DNase.t3", "H3K27ac.t3", "H3K4me1.t3",
              "H3K4me2.t3", "EP300.t3"]
TIME0_FEATURES=["DNase.t00", "H3K27ac.t00", "H3K4me1.t00", "H3K4me2.t00", 
              "EP300.t00"]
TIME3_FEATURES=["DNase.t3","H3K27ac.t3","H3K4me1.t3","H3K4me2.t3","EP300.t3"]
H3K2AC=["H3K27ac.t00", "H3K27ac.t3"]
P300=["EP300.t00", "EP300.t3"]
H3K4ME1=["H3K4me1.t00", "H3K4me1.t3"]
H3K4ME2=["H3K4me2.t00", "H3K4me2.t3"]
DNASE=["DNase.t00", "DNase.t3"]

def System(cmd):
    print(cmd)
    #os.system(cmd)

def getDTRK(features):
    keep=set()
    for feature in features: keep.add(feature)
    DTRK=""
    for feature in ALL_FEATURES:
        if(not feature in keep): DTRK+=" DTRK "+feature
    return DTRK

def subsetModel(originalModel,newModel,keepFeatures):
    DTRK=getDTRK(features)
    System("cp "+originalModel+" "+newModel)
    System(MUMMIE+"/hmm-edit "+newModel+" "+DTRK)

def retrainModel(model,fastb):
    TGF=DELTA+"/hmm/tgf.tgf"
    TIE=DELTA+"/hmm/tie.txt"
    temp="temp.hmm"
    cmd=BAUMWELCH\
        +" -S -W -t "+TIE+" -c 32 -L 0.01 -N 10000 "\
        +model+" "+TGF+" "+fastb+" 1000 "+temp
    System(cmd)
    System("mv "+temp+" "+model)

def subsetFastb(fastbDir,keepFeatures,outDir):
    cmd=MUMMIE+"/fastb-extract-tracks-dir.pl "+fastbDir+" "+outDir+" "\
        +keepFeatures
    System(cmd)

def getAUC(posHMM,negHMM,posFastb,negFastb)
    cmd=DELTA+"/src/delta-test.py "+posHMM+" "+negHMM+" "+posFASTB+" "\
        +negFASTB+" > roc.tmp; roc.pl roc.tmp > subset.roc"
    System(cmd)
    System("area-under-ROC.pl subset.roc > auc.txt")
    auc=loadAUC()

def loadAUC():
    with open("auc.txt","rt") as IN:
        line=IN.readline()
        auc=int(line.rstrip())
        return auc

def test(features,label,posHMM,negHMM,posTrain,negTrain,posTest,negTest):
    subsetModel(posHMM,TEMP_POS_HMM,features)
    subsetModel(negHMM,TEMP_NEG_HMM,features)
    retrainModel(TEMP_POS_HMM,posTrain)
    retrainModel(TEMP_NEG_HMM,negTrain)
    subsetFastb(posTest,features,TEMP_POS_TEST)
    subsetFastb(negTest,features,TEMP_NEG_TEST)
    getAUC(posHMM,negHMM,TEMP_POS_TEST,TEMP_NEG_TEST)

def testPartition(features,label,partition):
    i=partition-1
    posHMM=HMM_DIR+"/crossval-pos-bin"+partition+"-rep"+BEST_MODELS[i]+".hmm"
    negHMM=HMM_DIR+"/crossval-neg-bin"+partition+".hmm"
    posTrain=TRAIN_POS+"/leaveout"+partition
    negTrain=TRAIN_NEG+"/leaveout"+partition
    posTest=TEST_POS+"/"+partition
    negTest=TEST_NEG+"/"+partition
    test(ALL_FEATURES,"full_model",posHMM,negHMM,posTrain,negTrain,
         posTest,negTest)
    test(TIME0_FEATURES,"all_t00",posHMM,negHMM,posTrain,negTrain,
         posTest,negTest)
    test(TIME3_FEATURES,"all_t3",posHMM,negHMM,posTrain,negTrain,
         posTest,negTest)
    test(H3K2AC,"H3K2AC_both_times",posHMM,negHMM,posTrain,negTrain,
         posTest,negTest)
    test(P300,"P300_both_times",posHMM,negHMM,posTrain,negTrain,
         posTest,negTest)
    test(H3K4ME1,"H3K4ME1_both_times",posHMM,negHMM,posTrain,negTrain,
         posTest,negTest)
    test(H3K4ME2,"H3K4ME2_both_times",posHMM,negHMM,posTrain,negTrain,
         posTest,negTest)
    test(DNASE,"DNASE_both_times",posHMM,negHMM,posTrain,negTrain,
         posTest,negTest)
    test(["DNase.t00"],"DNase.t00",posHMM,negHMM,posTrain,negTrain,
         posTest,negTest)
    test(["H3K27ac.t00"],"H3K27ac.t00",posHMM,negHMM,posTrain,negTrain,
         posTest,negTest)
    test(["H3K4me1.t00"],"H3K4me1.t00",posHMM,negHMM,posTrain,negTrain,
         posTest,negTest)
    test(["H3K4me2.t00"],"H3K4me2.t00",posHMM,negHMM,posTrain,negTrain,
         posTest,negTest)
    test(["EP300.t00"],"EP300.t00",posHMM,negHMM,posTrain,negTrain,
         posTest,negTest)
    test(["DNase.t3"],"DNase.t3",posHMM,negHMM,posTrain,negTrain,
         posTest,negTest)
    test(["H3K27ac.t3"],"H3K27ac.t3",posHMM,negHMM,posTrain,negTrain,
         posTest,negTest)
    test(["H3K4me1.t3"],"H3K4me1.t3",posHMM,negHMM,posTrain,negTrain,
         posTest,negTest)
    test(["H3K4me2.t3"],"H3K4me2.t3",posHMM,negHMM,posTrain,negTrain,
         posTest,negTest)
    test(["EP300.t3"],"EP300.t3",posHMM,negHMM,posTrain,negTrain,
         posTest,negTest)

#=========================================================================
# main()
#=========================================================================

for i in range(NUM_PARTITIONS):
    partition=i+1
    testPartition(partition)


