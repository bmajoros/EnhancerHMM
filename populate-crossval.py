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
from SlurmWriter import SlurmWriter

GENERATE_PARTITIONS=False
MAKE_TRAINING_SETS=False
MAKE_SLURMS=True

NUM_PARTITIONS=5
REPLICATES=300
BASE="/home/bmajoros/GGR/delta/fastb"
POSITIVES=BASE+"/positives"
NEGATIVES=BASE+"/negatives"
CROSSVAL=BASE+"/crossval"
PARTITIONS=CROSSVAL+"/partitions"
PARTITIONS_POS=PARTITIONS+"/pos"
PARTITIONS_NEG=PARTITIONS+"/neg"
TRAINING_SETS=CROSSVAL+"/training-sets"
TRAINING_SETS_POS=TRAINING_SETS+"/pos"
TRAINING_SETS_NEG=TRAINING_SETS+"/neg"
MUMMIE=os.environ["MUMMIE"]
BAUMWELCH=MUMMIE+"/baum-welch"
HMM_DIR="/home/bmajoros/GGR/delta/hmm"
TEMPLATE_HMM=HMM_DIR+"/1path-best-full.hmm"
NEG_TEMPLATE_HMM=HMM_DIR+"/trained-neg-nomotif.hmm"
SLURM_DIR="/home/bmajoros/GGR/delta/slurms/crossval-slurms"

def populate(bin,partitions,trainingSets):
    dest=trainingSets+"/leaveout"+str(bin)
    for i in range(NUM_PARTITIONS):
        index=i+1
        if(index==bin): continue
        src=partitions+"/"+str(index)
        cmd="cp "+src+"/* "+dest
        os.system(cmd)

def partition(files):
    N=len(files)
    binSize=int(N/NUM_PARTITIONS)
    bins=[]
    for i in range(NUM_PARTITIONS): bins.append([])
    for i in range(N):
        bin=int(i/binSize)
        if(bin>=NUM_PARTITIONS): bin=NUM_PARTITIONS-1
        bins[bin].append(files[i])
    return bins

def copyBins(source,bins,dest):
    numBins=len(bins)
    for i in range(numBins):
        destDir=dest+"/"+str(i+1)
        print("COPY BIN "+str(i))
        copyBin(source,bins[i],destDir)

def copyBin(source,files,dest):
    for file in files:
        cmd="cp "+source+"/"+file+" "+dest
        os.system(cmd)

#=========================================================================
# main()
#=========================================================================

# First, generate partitions
if(GENERATE_PARTITIONS):
    positives=os.listdir(POSITIVES)
    negatives=os.listdir(NEGATIVES)
    posBins=partition(positives)
    negBins=partition(negatives)
    copyBins(POSITIVES,posBins,PARTITIONS_POS)
    copyBins(NEGATIVES,negBins,PARTITIONS_NEG)

# Merge partitions into training sets
if(MAKE_TRAINING_SETS):
    for i in range(NUM_PARTITIONS):
        populate(i+1,PARTITIONS_POS,TRAINING_SETS_POS)
        populate(i+1,PARTITIONS_NEG,TRAINING_SETS_NEG)

# Make slurm scripts for training
if(MAKE_SLURMS):
    writer=SlurmWriter()
    for i in range(NUM_PARTITIONS):
        bin=i+1
        traindir=TRAINING_SETS_NEG+"/leaveout"+str(bin)
        outfile=HMM_DIR+"/crossval-neg-bin"+str(bin)+".hmm"
        cmd="cd "+HMM_DIR+" ; "+BAUMWELCH+" -S -c 32 -L 0.01 -N 10000 "+NEG_TEMPLATE_HMM+" tgf.tgf "+traindir+" 1000 "+outfile
        writer.addCommand(cmd)
        for rep in range(REPLICATES):
            traindir=TRAINING_SETS_POS+"/leaveout"+str(bin)
            outfile=HMM_DIR+"/crossval-pos-bin"+str(bin)+"-rep"+str(rep)+".hmm"
            cmd="cd "+HMM_DIR+" ; "+BAUMWELCH+" -S -W -t tie.txt -c 32 -L 0.01 -N 10000 "+TEMPLATE_HMM+" tgf.tgf "+traindir+" 1000 "+outfile
            writer.addCommand(cmd)
writer.mem(5000)
writer.nice(500)
writer.threads(32)
writer.setQueue("new,all")
writer.writeArrayScript(SLURM_DIR,"CROSSVAL",SLURM_DIR,500)


