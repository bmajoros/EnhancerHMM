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
import ProgramName
from Rex import Rex
rex=Rex()

GENES_ENHANCERS="genes-enhancers-distal.txt"
CLUSTER_MEMBERS="figure/cluster-membership.txt"

#=========================================================================
# main()
#=========================================================================
clusterGeneCounts={}
clusterEnhancerCounts={}
enhancerToCluster={}
with open(CLUSTER_MEMBERS,"rt") as IN:
    for line in IN:
        fields=line.rstrip().split()
        if(len(fields)!=2): continue
        enhancer=fields[0]
        cluster=int(fields[1])
        if(not rex.find("\"(\S+)\"",enhancer)):
            raise Exception("no quotes")
        enhancer=rex[1]
        enhancerToCluster[enhancer]=cluster
        clusterEnhancerCounts[cluster]=clusterEnhancerCounts.get(cluster,0)+1
        clusterGeneCounts[cluster]=0
with open(GENES_ENHANCERS,"rt") as IN:
    for line in IN:
        fields=line.rstrip().split()
        if(len(fields)!=2): continue
        gene=fields[0]
        enhancer=fields[1]
        cluster=enhancerToCluster.get(enhancer)
        if(cluster is None): continue
        clusterGeneCounts[cluster]+=1
for cluster in clusterEnhancerCounts.keys():
    print(cluster,clusterEnhancerCounts[cluster],clusterGeneCounts[cluster])


