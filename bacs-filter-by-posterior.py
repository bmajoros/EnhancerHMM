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
import math
from BedReader import BedReader
from Pipe import Pipe
import SumLogProbs

MUMMIE=os.environ["MUMMIE"]

def getLL(fastb,hmm):
    cmd=MUMMIE+"/get-likelihood "+hmm+" "+fastb
    pipe=Pipe(cmd)
    line=pipe.readline()
    return float(line)

if(len(sys.argv)!=8):
    exit(sys.argv[0]+" <predictions.bed> <schema> <fg.hmm> <bg.hmm> <fg_prior> <bg_prior> <fastb-dir>")
(predictionsBed,schema,fgHMM,bgHMM,fgPrior,bgPrior,fastbDir)=sys.argv[1:]

fgPrior=math.log(float(fgPrior))
bgPrior=math.log(float(bgPrior))

reader=BedReader(predictionsBed)
while(True):
    record=reader.nextRecord() # Bed3Record or Bed6Record
    if(not record): break
    cmd=MUMMIE+"/fastb-slice "+schema+" "+fastbDir+"/"+record.chr+".t3.fastb "+str(record.interval.begin)+" "+str(record.interval.end)+" slice.fastb"
    os.system(cmd)
    fgLL=getLL("slice.fastb",fgHMM)
    bgLL=getLL("slice.fastb",bgHMM)
    jointBG=bgLL+bgPrior
    jointFG=fgLL+fgPrior
    denom=SumLogProbs.sumLogProbs2(jointBG,jointFG)
    posterior=math.exp(jointFG-denom)
    print(posterior)
reader.close()


