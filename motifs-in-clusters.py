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

if(len(sys.argv)!=4):
    exit(ProgramName.get()+" <hits.txt> <cluster-membership.txt> <cluster-sizes.txt>\n")
(hitsFile,membershipFile,clusterSizeFile)=sys.argv[1:]

clusterSizes={}
IN=open(clusterSizeFile,"rt")
i=0
for line in IN:
    clusterSizes[i]=int(line)
    i+=1
IN.close()

peakToCluster={}
IN=open(membershipFile,"rt")
for line in IN:
    if(not rex.find("\"(\S+)\"\t(\d+)",line)): continue
    peak=rex[1]
    cluster=int(rex[2])
    peakToCluster[peak]=cluster
IN.close()

counts=[0]*10
IN=open(hitsFile,"rt")
for line in IN:
    line=line.rstrip()
    if(not peakToCluster.get(line,None)): continue
    cluster=peakToCluster[line]-1
    counts[cluster]+=1
IN.close()

for i in range(10):
    size=clusterSizes[i]
    count=counts[i]
    proportion=float(count)/float(size)
    print(i+1,count,proportion,sep="\t")
