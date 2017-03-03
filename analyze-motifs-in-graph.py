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
from Rex import Rex
rex=Rex()

def regressSubset(graph):
    genes=graph.getGenes()
    for gene in genes:
        mates=gene.enhancerMates
        if(len(mates)<1): continue
        for enhancer in mates:
            #if(enhancer.numPeaks<2): continue
            if(enhancer.numPeaks>1): continue
            if(enhancerHasAP1(enhancer)): continue
            #if(not enhancerHasAP1(enhancer)): continue
            p300=getP300(enhancer)
            print(gene.logFC,p300,sep="\t")

def regressAll(graph):
    genes=graph.getGenes()
    for gene in genes:
        mates=gene.enhancerMates
        if(len(mates)<1): continue
        for enhancer in mates:
            #if(enhancer.numPeaks<2): continue
            if(enhancer.numPeaks>1): continue
            p300=getP300(enhancer)
            print(gene.logFC,p300,sep="\t")

def getP300(enhancer):
    p300=0
    for peak in enhancer.peaks:
        p300+=peak.raw_p300_t3-peak.raw_p300_t0
    return p300

def enhancerHasGR(enhancer):
    for peak in enhancer.peaks:
        if("GR" in peak.motifs): return True
    return False

def enhancerHasAP1(enhancer):
    for peak in enhancer.peaks:
        for motif in peak.motifs:
            if("JUN" in motif or "jun" in motif or
               "FOS" in motif or "fos" in motif): return True
    return False

def motifsInNonGRpeaks(graph):
    enhancers=graph.getEnhancers()
    TFs=graph.getAllTFs()
    for tf in TFs:
        table=[[0,0],[0,0]] # [motif present][GR present)
        for enhancer in enhancers:
            #if(enhancer.numPeaks<2): continue
            if(enhancer.numPeaks>1): continue
            for peak in enhancer.peaks:
                x=1 if "GR" in peak.motifs else 0
                #found=0
                #for motif in peak.motifs:
                #    if("JUN" in motif or "jun" in motif or
                #       "FOS" in motif or "fos" in motif): 
                #        found=True
                #        break
                #x=1 if found else 0
                y=1 if tf in peak.motifs else 0
                table[x][y]+=1
        fisher=Fisher2x2(table[0][0],table[0][1],table[1][0],table[1][1])
        P=fisher.getPvalue()
        (exp00,exp01,exp10,exp11)=fisher.getExpectedCounts()
        if(P>=0.05): continue
        #if(table[0][1]<=exp01):continue
        msg="with GR" if table[1][1]>exp11 else "non-GR"
        print(P,tf,"observed=["+str(table[0][0])+" "+str(table[0][1])+" "+
              str(table[1][0])+" "+str(table[1][1])+"] expected=["+str(exp00)+
              " "+str(exp01)+" "+str(exp10)+" "+str(exp11)+"]",msg,
              sep="\t",flush=True)

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
        #if(P>=0.05): continue
        #if(table[1][0]<=exp10):continue
        #P=table[1][0]-exp10
        P=table[1][0]
        msg="DEX responsive" if table[1][1]>exp11 else "non-reponsive"
        print(P,tf,"observed=["+str(table[0][0])+" "+str(table[0][1])+" "+
              str(table[1][0])+" "+str(table[1][1])+"] expected=["+str(exp00)+
              " "+str(exp01)+" "+str(exp10)+" "+str(exp11)+"]",msg,
              sep="\t",flush=True)

def GR_vs_AP1(graph):
    GR_AP1(graph.getSingletons(),"singletons")
    GR_AP1(graph.getMultipeakEnhancers(),"multipeaks")

def GR_AP1(enhancers,label):
    gr=round(getProportionWithMotif(enhancers,"GR"),3)
    ap1=round(getProportionWithAP1(enhancers),3)
    both=round(getProportionWithGR_AP1(enhancers),3)
    print(label,"GR="+str(gr),"AP-1="+str(ap1),"both="+str(both),sep="\t")

def getProportionWithMotif(enhancers,motif):
    total=0
    hits=0
    for enhancer in enhancers:
        total+=1
        found=False
        for peak in enhancer.peaks:
            #total+=1
            #found=False
            if(motif in peak.motifs):
                found=True
            #if(found): hits+=1
        if(found): hits+=1
    proportion=float(hits)/float(total)
    return proportion

def getProportionWithAP1(enhancers):
    total=0
    hits=0
    for enhancer in enhancers:
        total+=1
        found=False
        for peak in enhancer.peaks:
            #total+=1
            #found=False
            for motif in peak.motifs:
                if(rex.find("JUN",motif) or rex.find("jun",motif) or \
                       rex.find("FOS",motif) or rex.find("fos",motif)):
                    found=True
            #if(found): hits+=1
        if(found): hits+=1
    proportion=float(hits)/float(total)
    return proportion

def getProportionWithGR_AP1(enhancers):
    total=0
    hits=0
    for enhancer in enhancers:
        total+=1
        foundGR=False; foundAP1=False
        for peak in enhancer.peaks:
            #total+=1
            #foundGR=False; foundAP1=False
            if("GR" in peak.motifs): foundGR=True
            for motif in peak.motifs:
                if(rex.find("JUN",motif) or rex.find("jun",motif) or
                   rex.find("FOS",motif) or rex.find("fos",motif)):
                    foundAP1=True
            #if(foundGR and foundAP1): hits+=1
        if(foundGR and foundAP1): hits+=1
    proportion=float(hits)/float(total)
    return proportion

#=========================================================================
# main()
#=========================================================================
if(len(sys.argv)!=3):
    exit(ProgramName.get()+" <enhancer-timepoint> <gene-timepoint>\n")
(enhancerTimepoint,geneTimepoint)=sys.argv[1:]

graph=EnhancerGraph(enhancerTimepoint,geneTimepoint)
#multiResponsiveMotifs(graph)
#GR_vs_AP1(graph)
#motifsInNonGRpeaks(graph)
regressSubset(graph)

