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

GENES_ENHANCERS="/home/bmajoros/GGR/delta/genes-enhancers-distal.txt"
#  ENSG00000176055.9       iter0_peak14463
CLUSTER_FILE="/home/bmajoros/GGR/delta/figure/cluster-membership.txt"
#  "iter0_peak843" 10

def saveCluster(genes,clusterID):
    filename="/home/bmajoros/GGR/delta/gene-cluster-"+str(clusterID)+".txt"
    with open(filename,"wt") as OUT:
        for geneID in genes:
            print(geneID,file=OUT)

#=========================================================================
# main()
#=========================================================================

# Load cluster membership file for enhancers
enhancerToCluster={}
largestClusterID=0
with open(CLUSTER_FILE,"rt") as IN:
    for line in IN:
        fields=line.rstrip().split()
        if(len(fields)!=2): continue
        (enhancerID,clusterID)=fields
        if(not rex.find("\"(\S+)\"",enhancerID)):
            raise Exception("Can't parse: "+enhancerID)
        enhancerID=rex[1]
        clusterID=int(clusterID)
        enhancerToCluster[enhancerID]=clusterID
        if(clusterID>largestClusterID): largestClusterID=clusterID
        
# Load the pairing between genes and enhancers, to put genes into clusters
genesByCluster=[]
for i in range(largestClusterID+1): genesByCluster.append([])
with open(GENES_ENHANCERS,"rt") as IN:
    for line in IN:
        fields=line.rstrip().split()
        if(len(fields)!=2): continue
        (geneID,enhancerID)=fields
        if(rex.find("(\S+)\.\d+",geneID)): geneID=rex[1]
        clusterID=enhancerToCluster.get(enhancerID,None)
        if(clusterID is None):
            continue
            #raise Exception("No cluster for "+enhancerID)
        genesByCluster[clusterID].append(geneID)
        
# Print genes by cluster
for i in range(1,largestClusterID+1):
    genes=genesByCluster[i]
    saveCluster(genes,i)
    #print(i,end="")
    #for geneID in genes: print("\t",geneID,end="")
    #print()

