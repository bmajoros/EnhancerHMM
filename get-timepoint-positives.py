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
import ProgramName
from Rex import Rex
rex=Rex()

#=========================================================================
# main()
#=========================================================================
if(len(sys.argv)!=5):
    exit(ProgramName.get()+" <genomewide-1path.txt> <delta-fastb> <min-LLR> <outdir-base>\n")
(predictionsFile,indir,minLLR,outbase)=sys.argv[1:]
MIN_LLR=float(minLLR)

with open(predictionsFile) as IN:
    for line in IN:
        fields=line.split()
        (fastb,begin,end,taskID,LLR,posterior)=fields[:6]
        if(float(LLR)<MIN_LLR): continue
        if(not rex.find("(\S+)\.(t\d+)\.fastb",fastb)):
            raise Exception("can't parse fastb: "+fastb)
        base=rex[1]; time=rex[2]
        outfile=outbase+"/"+time+"/"+fastb
        cmd="cp "+indir+"/"+fastb+" "+outfile
        #print(cmd)
        os.system(cmd)

