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
from Pipe import Pipe
import SumLogProbs
import os
import math

BASE="/home/bmajoros/GGR/p300"
SLICED=BASE+"/sliced"
PEAKS=BASE+"/peaks.txt"
POS_MODEL=BASE+"/model/sliced.hmm"
NEG_MODEL=BASE+"/model/trained-neg1.hmm"
MUMMIE=os.environ["MUMMIE"]

def getLikelihood(filename,hmm):
    cmd=MUMMIE+"/get-likelihood "+hmm+" "+filename
    pipe=Pipe(cmd)
    line=pipe.readline()
    line=line.rstrip()
    LL=float(line)
    return LL

def getPosterior(filename):
    posLL=getLikelihood(filename,POS_MODEL)
    negLL=getLikelihood(filename,NEG_MODEL)
    uniformPrior=math.log(0.5)
    Z=SumLogProbs.sumLogProbs2(posLL+uniformPrior,negLL+uniformPrior)
    posterior=posLL+uniformPrior-Z
    return math.exp(posterior)

def processPeak(peak):
    file0=SLICED+"/"+peak+".standardized_across_all_timepoints.t00.fastb"
    file5=SLICED+"/"+peak+".standardized_across_all_timepoints.t05.fastb"
    post0=getPosterior(file0)
    post5=getPosterior(file5)
    post0=int(post0*100.0+5.0/9.0)/100.0
    post5=int(post5*100.0+5.0/9.0)/100.0
    print(peak,post0,post5,sep="\t",flush=True)

#=============================== main() ===============================
with open(PEAKS,"rt") as fh:
    while(True):
        line=fh.readline()
        if(not line): break
        line=line.rstrip()
        processPeak(line)

    

