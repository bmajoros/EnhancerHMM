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
import ProgramName

TIMES=(0.5,1,2,3,4,5,6,7,8,10,12)

def getProportions(trajectories,count):
    prop=[0]*11
    #N=float(len(trajectories))
    N=[0]*11
    for traj in trajectories:
        for i in range(len(traj)):
            if(int(traj[i])>0): N[i]+=1
    for traj in trajectories:
        for i in range(len(traj)):
            if(int(traj[i])==count):
                #prop[i]+=1.0/N[i]
                prop[i]+=1.0
    return prop

def write(trajectory,filename):
    with open(filename,"wt") as OUT:
        for i in range(len(trajectory)):
            OUT.write(str(TIMES[i])+"\t"+str(trajectory[i])+"\n")

#=========================================================================
# main()
#=========================================================================
if(len(sys.argv)!=4):
    exit(ProgramName.get()+" <trajectories.txt> <out-single.txt> <out-dual.txt>\n")
(trajectoryFile,singleFile,dualFile)=sys.argv[1:]

trajectories=[]
with open(trajectoryFile,"rt") as IN:
    for line in IN:
        fields=line.rstrip().split()
        if(len(fields)!=12): continue
        trajectory=fields[1:]
        trajectories.append(trajectory)
single=getProportions(trajectories,1)
dual=getProportions(trajectories,2)
write(single,singleFile)
write(dual,dualFile)

