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
import sys
import ProgramName
from Rex import Rex
rex=Rex()

SUBSET_ENHANCERS=set()
SHOULD_SUBSET=True

def loadLoopFile(filename):
    hash={} # maps genes to enhancers
    with open(filename) as IN:
        for line in IN:
            fields=line.rstrip().split()
            if(len(fields)!=2): continue
            (gene,peak)=fields
            if(hash.get(gene,None) is None): hash[gene]=[]
            hash[gene].append(peak)
    return hash

def loadExpressionFile(filename):
    # gene        logFC   logCPM  LR      PValue  FDR
    hash={}
    with open(filename) as IN:
        for line in IN:
            fields=line.rstrip().split()
            if(len(fields)!=6): continue
            (gene,logFC,logCPM,LR,Pvalue,FDR)=fields
            hash[gene]=logFC
    return hash

#=========================================================================
# main()
#=========================================================================
if(len(sys.argv)!=5):
    exit(ProgramName.get()+" <genes-enhancers.txt> <edgeR.txt> <cluster-membership.txt> <cluster-number>\n")
(loopFile,expressionFile,membershipFile,whichCluster)=sys.argv[1:]
whichCluster=int(whichCluster)

if(not os.path.exists(membershipFile)): exit(membershipFile+" not found")
with open(membershipFile) as IN:
    for line in IN:
        fields=line.rstrip().split()
        if(len(fields)!=2): continue
        (id,cluster)=fields
        if(int(cluster)!=whichCluster): continue
        if(not rex.find("\"(\S+)\"",id)): exit(id)
        SUBSET_ENHANCERS.add(rex[1])
        
geneToEnhancer=loadLoopFile(loopFile)
geneToFC=loadExpressionFile(expressionFile)

genes=geneToFC.keys()
for gene in genes:
    enhancers=geneToEnhancer.get(gene,None)
    if(enhancers is None or len(enhancers)<1): continue
    found=False
    for enhancer in enhancers:
        #print("checking",gene,enhancer)
        if(enhancer in SUBSET_ENHANCERS): found=True
    if(not found): continue
    Y=geneToFC[gene]
    print(Y)



