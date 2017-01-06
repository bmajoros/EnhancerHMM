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
import TempFilename
from Fastb import Fastb

BASE="/home/bmajoros/GGR/delta"
HMM_DIR=BASE+"/hmm"
POS_HMM=HMM_DIR+"/trained-pos5.hmm"
NEG_HMM=HMM_DIR+"/trained-neg1.hmm"
TEST_POS=BASE+"/test-pos"
TEST_NEG=BASE+"/test-neg"
MUMMIE=os.environ["MUMMIE"]
EVAL=MUMMIE+"/get-likelihood"
TEMP=TempFilename.generate(".fastb")

def evaluate(dir,label):
    files=os.listdir(dir)
    for file in files:
        path=dir+"/"+file
        fastb=Fastb(path)
        track=fastb.getTrackByName("EP300.t3")
        (maxValue,maxPos)=track.getMax()
        print(maxValue,label,sep="\t")

evaluate(TEST_POS,1)
evaluate(TEST_NEG,0)
os.remove(TEMP)


