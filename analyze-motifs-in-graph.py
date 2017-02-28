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
from EnhancerGraph import EnhancerGraph
from Fisher2x2 import Fisher2x2

def multiResponsiveMotifs(graph):
    enhancers=graph.getEnhancers()
    TFs=graph.getAllTFs()
    for tf in TFs:
        table=[[0,0],[0,0]] # [motif present][dex responsive]
        for enhancer in enhancers:
            if(enhancer.numPeaks<2): continue
            for peak in enhancer.peaks:
                x=1 if tf in peak.motifs else 0
                y=1 if peak.dex_responsive else 0
                table[x][y]+=1
        fisher=Fisher2x2(table[0][0],table[0][1],table[1][0],table[1][1])
        P=fisher.getPvalue()
        (exp00,exp01,exp10,exp11)=fisher.getExpectedCounts()
        if(P>=0.05): continue
        msg="DEX responsive" if table[1][1]>exp11 else "non-reponsive"
        print(P,tf,"observed=["+str(table[0][0])+" "+str(table[0][1])+" "+
              str(table[1][0])+" "+str(table[1][1])+"] expected=["+str(exp00)+
              " "+str(exp01)+" "+str(exp10)+" "+str(exp11)+"]",msg,
              sep="\t",flush=True)

#=========================================================================
# main()
#=========================================================================
if(len(sys.argv)!=3):
    exit(ProgramName.get()+" <enhancer-timepoint> <gene-timepoint>\n")
(enhancerTimepoint,geneTimepoint)=sys.argv[1:]

graph=EnhancerGraph(enhancerTimepoint,geneTimepoint)
multiResponsiveMotifs(graph)



