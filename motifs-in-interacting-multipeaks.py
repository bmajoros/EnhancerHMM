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
from Fisher2x2 import Fisher2x2
from EssexParser import EssexParser

MAX_ELEMENTS=1000000000

def buildGraph(infile):
    genes=[]; enhancers=[]
    parser=EssexParser(infile)
    numRead=0
    while(True):
        node=parser.nextElem()
        if(node is None): break
        if(node.getTag()=="gene"): genes.append(node)
        elif(node.getTag()=="enhancer"): enhancers.append(node)
        else: raise Exception(node.getTag())
        numRead+=1
        if(numRead>=MAX_ELEMENTS): break
    linkObjects(genes,enhancers)
    return (genes,enhancers)

def linkObjects(genes,enhancers):
    geneHash={}; enhancerHash={}
    hashIDs(genes,geneHash)
    hashIDs(enhancers,enhancerHash)
    linkup(genes,geneHash,enhancerHash)
    linkup(enhancers,geneHash,enhancerHash)

def hashIDs(array,hash):
    for elem in array: hash[elem.getAttribute("id")]=elem

def linkup(array,genes,enhancers):
    for elem in array:
        elem.geneMates=[]; elem.enhancerMates=[]
        matesNode=elem.findChild("mates")
        for mateID in matesNode.elements:
            if(genes.get(mateID,None) is not None):
                elem.geneMates.append(genes[mateID])
            elif(enhancers.get(mateID,None) is not None):
                elem.enhancerMates.append(enhancers[mateID])

def getConnectedMultis(enhancers):
    keep=[]
    for enhancer in enhancers:
        if(int(enhancer.getAttribute("numPeaks"))<2): continue
        if(len(enhancer.enhancerMates)>0):
            enhancer.connected=True
            enhancer.motifs=set()
            keep.append(enhancer)
    return keep

def getDisconnectedMultis(enhancers):
    keep=[]
    for enhancer in enhancers:
        if(int(enhancer.getAttribute("numPeaks"))<2): continue
        if(len(enhancer.enhancerMates)==0):
            enhancer.connected=False
            enhancer.motifs=set()
            keep.append(enhancer)
    return keep

def getAllTFs(enhancers):
    TFs=set()
    for enhancer in enhancers:
        motifNodes=enhancer.findDescendents("motifs")
        for motifNode in motifNodes:
            for motif in motifNode.elements:
                TFs.add(motif)
                enhancer.motifs.add(motif)
                #enhancer.motifs.add("ANY_MOTIF")
    return TFs

def testTF(tf,allMultis):
    table=[[0,0],[0,0]]
    for enhancer in allMultis:
        x=1 if tf in enhancer.motifs else 0
        y=1 if enhancer.connected else 0
        table[x][y]+=1
    fisher=Fisher2x2(table[0][0],table[0][1],table[1][0],table[1][1])
    P=fisher.getPvalue()
    if(P>0.05): return
    (exp00,exp01,exp10,exp11)=fisher.getExpectedCounts()
    msg="enriched in connected" if table[1][1]>exp11 \
        else "enriched in disconnected"
    print(P,tf,"observed=["+
          str(table[0][0])+" "+str(table[0][1])+" "+
          str(table[1][0])+" "+str(table[1][1])+"] expected=["+
          str(exp00)+" "+str(exp01)+" "+str(exp10)+" "+str(exp11)+"]",
          msg,sep="\t")

#=========================================================================
# main()
#=========================================================================
if(len(sys.argv)!=2):
    exit(ProgramName.get()+" <in.essex>\n")
(infile,)=sys.argv[1:]

(genes,enhancers)=buildGraph(infile)
connected=getConnectedMultis(enhancers)
disconnected=getDisconnectedMultis(enhancers)
allMultis=[]; allMultis.extend(connected); allMultis.extend(disconnected)
TFs=getAllTFs(allMultis)
for TF in TFs: testTF(TF,allMultis)
