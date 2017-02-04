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
import os
from Pipe import Pipe

BASE="/home/bmajoros/GGR/delta"
#HMM_DIR=BASE+"/hmm"
#POS_HMM=HMM_DIR+"/1path-best-half.hmm"
#NEG_HMM=HMM_DIR+"/trained-neg1.hmm"
#POS_HMM=HMM_DIR+"/subset-pos.hmm"
#NEG_HMM=HMM_DIR+"/subset-neg.hmm"
#POS_HMM=HMM_DIR+"/trained-twopaths.hmm"
#POS_HMM=HMM_DIR+"/twopaths-bestLL.hmm"
#NEG_HMM=HMM_DIR+"/trained-bg-twopaths.hmm"
#POS_HMM=HMM_DIR+"/onepath.hmm"
#POS_HMM=HMM_DIR+"/retrain.hmm"
#NEG_HMM=HMM_DIR+"/trained-neg-nomotif.hmm"
#POS_HMM=HMM_DIR+"/2path-62.hmm"
#POS_HMM=HMM_DIR+"/twopaths-dna2x.hmm"
#NEG_HMM=HMM_DIR+"/neg-dna2x.hmm"
#TEST_POS=BASE+"/subset-pos"
#TEST_NEG=BASE+"/subset-neg"
#TEST_POS=BASE+"/test-pos-nomotif"
#TEST_NEG=BASE+"/test-neg-nomotif"
#TEST_POS=BASE+"/test-pos-dna"
#TEST_NEG=BASE+"/test-neg-dna"
MUMMIE=os.environ["MUMMIE"]
EVAL=MUMMIE+"/get-likelihood"

def getLL(file,hmm):
    cmd=EVAL+" "+hmm+" "+file
    pipe=Pipe(cmd)
    line=pipe.readline()
    #print(cmd)
    LL=float(line.rstrip())
    return LL

def evaluate(dir,label):
    files=os.listdir(dir)
    for file in files:
        path=dir+"/"+file
        posLL=getLL(path,POS_HMM)
        negLL=getLL(path,NEG_HMM)
        LLR=posLL-negLL
        print(LLR,label,sep="\t")

#=========================================================================
# main()
#=========================================================================
if(len(sys.argv)!=5):
    exit(ProgramName.get()+" <pos.hmm> <neg.hmm> <test-pos> <test-neg>\n")
(POS_HMM,NEG_HMM,TEST_POS,TEST_NEG)=sys.argv[1:]

evaluate(TEST_POS,1)
evaluate(TEST_NEG,0)



