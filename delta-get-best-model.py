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
import sys
import ProgramName
from Rex import Rex
rex=Rex()

#=========================================================================
# main()
#=========================================================================
if(len(sys.argv)!=3):
    exit(ProgramName.get()+" <outputs-dir> <top-N>\n")
(BASE,N)=sys.argv[1:]

#BASE="/home/bmajoros/GGR/delta/slurms/train-slurms/outputs"
#BASE="/home/bmajoros/GGR/delta/slurms/train-full-slurms/outputs"
files=os.listdir(BASE)
models=[]
for file in files:
    LL=None
    done=False
    with open(BASE+"/"+file,"rt") as IN:
        for line in IN:
            if(rex.find("ITERATION.*LL=(\S+)\s+deltaLL",line)):
                LL=float(rex[1])
            elif(rex.find("^done",line)): done=True
    if(not done): continue
    if(LL is not None):
        if(not rex.find("(\d+).output",file)): exit("can't parse: "+file)
        index=int(rex[1])
        bestModel=str(index)+".hmm"
        models.append([bestModel,LL])
models.sort(key=lambda x: -x[1])
for i in range(int(N)):
    pair=models[i]
    (file,LL)=pair
    print(file,LL,sep="\t")


