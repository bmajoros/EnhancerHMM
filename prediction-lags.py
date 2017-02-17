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
from Rex import Rex
rex=Rex()

TIMEPOINTS=("t05","t1","t2","t3","t4","t5","t6","t7","t8","t10","t12")

for time in TIMEPOINTS: print("\t"+time,end="")
print()
N=len(TIMEPOINTS)
for i in range(N):
#for i in range(N-3,N):
    enhancerTime=TIMEPOINTS[i]
    print(enhancerTime+"\t",end="")
    for j in range(N):
        geneTime=TIMEPOINTS[j]
        cmd="src/regress-on-total-P300.py genes-enhancers-distal.txt genomewide-features-raw.txt /data/reddylab/projects/GGR/results/rna_seq/differential_expression/iter0/edgeR/edgeR.sva."+geneTime+".vs.t00.protein_coding.txt "+enhancerTime+" > regress-total-p300.txt"
        os.system(cmd)
        pipe=Pipe("lm.R regress-total-p300.txt")
        while(True):
            line=pipe.readline()
            if(line is None): break
            if(rex.find("Adjusted R-squared:\s*(\S+)",line)):
                R=float(rex[1])
                R=round(R,2)
                print(str(R)+"\t",end="")
                break
    print()




