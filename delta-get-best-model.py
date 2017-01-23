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
import os
from Rex import Rex
rex=Rex()

BASE="/home/bmajoros/GGR/delta/slurms/train-slurms/outputs"
files=os.listdir(BASE)
bestModel=None
bestLL=None
for file in files:
    LL=None
    done=False
    with open(BASE+"/"+file,"rt") as IN:
        for line in IN:
            if(rex.find("ITERATION.*LL=(\S+)\s+deltaLL",line)):
                LL=float(rex[1])
            elif(rex.find("^done",line)): done=True
    #if(not done): continue
    if(bestLL is None or LL is not None and LL>bestLL):
        bestLL=LL
        if(not rex.find("(\d+).output",file)):
            exit("can't parse filename: "+file)
        index=int(rex[1])
        bestModel="threepaths-trained"+str(index-1)+".hmm"
print(bestModel,bestLL)


