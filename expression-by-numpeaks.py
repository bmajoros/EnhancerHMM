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

def getExpressions(enhancers,filename):
    with open(filename,"wt") as OUT:
        for enhancer in enhancers:
            for mate in enhancer.mates:
                if(mate.type!="gene"): continue
                print(mate.logFC,file=OUT)

def writeExpressions(genes,filename):
    with open(filename,"wt") as OUT:
        for gene in genes:
            print(gene.logFC,file=OUT)

def genesWithSingleton(singletons):
    genes=set()
    for enhancer in singletons:
        for mate in enhancer.mates:
            if(mate.type=="gene"): genes.add(mate)
    return genes

def genesWithNoSingleton(multipeaks,genesWithSingletons):
    genes=set()
    for enhancer in multipeaks:
        for mate in enhancer.mates:
            if(mate.type=="gene" and mate not in genesWithSingletons):
                genes.add(mate)
    return genes
            

#=========================================================================
# main()
#=========================================================================
if(len(sys.argv)!=3):
    exit(ProgramName.get()+" <enhancer-timepiont> <gene-timepoint>\n")
(enhancerTimepoint,geneTimepoint)=sys.argv[1:]

graph=EnhancerGraph(enhancerTimepoint,geneTimepoint)
#   genes=graph.getGenes()
#   enhancers=graph.getEnhancers()
singletons=graph.getSingletons()
multipeakEnhancers=graph.getMultipeakEnhancers()
#getExpressions(singletons,"singletons.tmp")
#getExpressions(multipeakEnhancers,"multipeaks.tmp")
genes1=genesWithSingleton(singletons)
genes2=genesWithNoSingleton(multipeakEnhancers,genes1)
writeExpressions(genes1,"singletons.tmp")
writeExpressions(genes2,"multipeaks.tmp")

