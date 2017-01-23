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
import ProgramName
from Rex import Rex
rex=Rex()

def readComponent(IN):
    means=IN.readline().rstrip().split()[1:]
    means=list(map(float,means))
    matrixSize=int(IN.readline().rstrip().split()[0])
    cov=[]
    for i in range(matrixSize):
        cov.append(IN.readline().rstrip())
    IN.readline()
    return [means,cov]

def whichComponent(weights):
    for i in range(len(weights)):
        if(float(weights[i])>0): return i
    return None

def parseState(state,IN):
    numComponents=int(IN.readline())
    weights=IN.readline().rstrip().split()
    components=[]
    for i in range(numComponents):
        component=readComponent(IN)
        components.append(component)
    index=whichComponent(weights)
    return components[index]

def emit0and3hr(parms,label):
    emitFeatures(range(0,5),parms,label+"_0hr",0)
    emitFeatures(range(5,10),parms,label+"_3hr",3)

def emitDiff(parms,label):
    emitDiffFeatures(range(0,5),parms,label)

def emitDiffFeatures(features,parms,label):
    outfile=label+".diffs"
    OUT=open(outfile,"wt")
    for f in features:
        for i in range(len(parms)):
            (state,component)=parms[i]
            (means,cov)=component
            print(means[f+5]-means[f],end="",file=OUT)
            print("\t" if i+1<len(parms) else "\n",end="",file=OUT)
    OUT.close()

def emitFeatures(features,parms,label,hour):
    outfile=label+".means"
    OUT=open(outfile,"wt")
    for f in features:
        for i in range(len(parms)):
            (state,component)=parms[i]
            (means,cov)=component
            print(means[f],end="",file=OUT)
            print("\t" if i+1<len(parms) else "\n",end="",file=OUT)
            covFile="state"+str(state)+"."+str(hour)+"hr.cov"
            writeCovFile(cov,covFile)
    OUT.close()

def writeCovFile(cov,filename):
    OUT=open(filename,"wt")
    for line in cov:
        print(line,file=OUT)
    OUT.close()

#=========================================================================
# main()
#=========================================================================
if(len(sys.argv)!=2):
    exit(ProgramName.get()+" <in.hmm>\n")
hmmFile=sys.argv[1]

# Read the HMM file
IN=open(hmmFile,"rt")
state=0
parms=[]
for line in IN:
    if(rex.find("state (\d+) emissions:",line)):
        state=rex[1]
        component=parseState(state,IN)
        parms.append([state,component])
IN.close()

# Generate output
numStates=len(parms)
if(numStates==10):
    emit0and3hr(parms[0:5],"path1")
    emit0and3hr(parms[5:10],"path2")
    emitDiff(parms[0:5],"path1")
    emitDiff(parms[5:10],"path2")
if(numStates==15):
    emit0and3hr(parms[0:5],"path1")
    emit0and3hr(parms[5:10],"path2")
    emit0and3hr(parms[10:15],"path3")

