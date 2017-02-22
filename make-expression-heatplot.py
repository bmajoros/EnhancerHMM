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
import os
from Pipe import Pipe
from SummaryStats import SummaryStats
from Rex import Rex
rex=Rex()

TIMEPOINTS=("t05","t1","t2","t3","t4","t5","t6","t7","t8","t10","t12")

N=len(TIMEPOINTS)
for cluster in range(1,11):
    first=True
    for geneTime in TIMEPOINTS:
        cmd="src/cluster-to-gene-expression.py genes-enhancers-distal.txt /data/reddylab/projects/GGR/results/rna_seq/differential_expression/iter0/edgeR/edgeR.sva."+geneTime+".vs.t00.protein_coding.txt "+" figure/cluster-membership.txt "+str(cluster)
        pipe=Pipe(cmd)
        array=[]
        while(True):
            line=pipe.readline()
            if(line is None): break
            if(not rex.find("\d",line)): continue
            array.append(float(line))
        (mean,SD,min,max)=SummaryStats.roundedSummaryStats(array)
        if(not first): print("\t",end="")
        print(mean,end="")
        first=False
    print()




