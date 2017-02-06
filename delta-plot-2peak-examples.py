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

#TIMEPOINTS=("t00","t05","t1","t2")
TIMEPOINTS=("t05",)
MUMMIE=os.environ["MUMMIE"]

def emit(id,fastbDir,timepoint):
    #filename=fastbDir+"/"+id+".standardized_across_all_timepoints."+timepoint+\
    #    ".fastb"
    filename=fastbDir+"/"+id+"."+timepoint+".fastb"
    cmd="cat "+filename+" | "+MUMMIE+"/fastb-to-xgraph.pl -t "+timepoint
    print(filename)
    os.system(cmd)
    

#=========================================================================
# main()
#=========================================================================
if(len(sys.argv)!=3):
    exit(ProgramName.get()+" <trajectories-num-peaks.txt> <fastb-dir>\n")
(trajectoryFile,fastb)=sys.argv[1:]

IN=open(trajectoryFile,"rt")
for line in IN:
    fields=line.rstrip().split()
    if(len(fields)!=12): continue
    (id,halfHour,oneHour)=fields[:3]
    halfHour=int(halfHour); oneHour=int(oneHour)
    #if(halfHour!=2): continue
    if(halfHour!=1): continue
    for t in TIMEPOINTS:
        emit(id,fastb,t)
IN.close()


