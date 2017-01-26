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
from Rex import Rex
rex=Rex()

# Process command line
if(len(sys.argv)!=4):
    exit(ProgramName.get()+" <cluster-membership.txt> <trajectories.txt> <1|2>\n")
(clusterFile,trajectoryFile,NUM_PEAKS)=sys.argv[1:]
NUM_PEAKS=int(NUM_PEAKS)

# Read trajectories
trajectories={}
IN=open(trajectoryFile,"rt")
for line in IN:
    fields=line.rstrip().split()
    if(len(fields)!=12): continue
    peak=fields[0]
    trajectories[peak]=fields[1:]
IN.close()

# Read cluster membership
clusters=[]
for i in range(10): clusters.append([])
IN=open(clusterFile,"rt")
for line in IN:
    fields=line.rstrip().split()
    if(len(fields)!=2): continue
    (peak,cluster)=fields
    if(not rex.find("\"(\S+)\"",peak)): exit("can't parse "+peak)
    peak=rex[1]
    cluster=int(cluster)-1
    clusters[cluster].append(peak)
IN.close()

# For each cluster, compute its centroid
for cluster in clusters:
    centroid=[0]*11
    N=0
    for peak in cluster:
        N+=1
        traj=trajectories[peak]
        for i in range(len(traj)):
            if(int(traj[i])==NUM_PEAKS): centroid[i]+=1
            #print(int(traj[i]),"versus",NUM_PEAKS)
    #print("N=",N)
    #print(centroid)
    for i in range(11): centroid[i]/=float(N)
    for i in range(11):
        print(centroid[i],end="\t" if i+1<11 else "\n")

