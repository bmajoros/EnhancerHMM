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
from Interval import Interval
from EssexNode import EssexNode
from Rex import Rex
rex=Rex()

CHUNK_WIDTH=1000
MAX_DIST_TO_ANCHOR=1

def loadLoops(filename):
    byChr={}
    IN=open(filename,"rt")
    for line in IN:
        fields=line.rstrip().split()
        if(len(fields)!=6): continue
        (fromChr,fromBegin,fromEnd,toChr,toBegin,toEnd)=fields
        fromBegin=int(fromBegin); fromEnd=int(fromEnd)
        toBegin=int(toBegin); toEnd=int(toEnd)
        interval1=makeLoopAnchor(fromChr,fromBegin,fromEnd)
        interval2=makeLoopAnchor(toChr,toBegin,toEnd)
        interval1.mate=interval2
        interval2.mate=interval1
        addAnchor(interval1,byChr)
        addAnchor(interval2,byChr)
    IN.close()
    return byChr

def addAnchor(interval,hash):
    if(hash.get(interval.chr,None) is None): hash[interval.chr]=[]
    hash[interval.chr].append(interval)

def makeLoopAnchor(fromChr,fromBegin,fromEnd):
    interval=Interval(fromBegin,fromEnd)
    interval.chr=fromChr
    interval.type="anchor"
    interval.objects=set()
    return interval

def loadExpressionFile(filename):
    hash={}
    with open(filename) as IN:
        for line in IN:
            fields=line.rstrip().split()
            if(len(fields)!=6): continue
            (gene,logFC,logCPM,LR,Pvalue,FDR)=fields
            hash[gene]=[logFC,logCPM,LR,Pvalue,FDR]
    return hash

def loadTssFile(filename,byChr,expressionHash):
    IN=open(filename,"rt")
    for line in IN:
        fields=line.rstrip().split()
        if(len(fields)!=4): continue
        (chr,pos,strand,id)=fields
        if(chr=="chrM"): continue
        pos=int(pos)
        interval=Interval(pos,pos+1)
        interval.chr=chr
        interval.type="gene"
        interval.id=id
        interval.mates=[]
        if(expressionHash.get(id,None) is None): continue
        (logFC,logCPM,LR,Pvalue,FDR)=expressionHash[id]
        interval.logFC=round(float(logFC),3)
        interval.FDR=round(float(FDR),3)
        byChr[chr].append(interval)
    IN.close()

def loadEnhancers(filename,byChr,byID):
    IN=open(filename,"rt")
    hash={}
    for line in IN:
        fields=line.rstrip().split()
        if(len(fields)!=4): continue
        (chr,begin,end,id)=fields
        if(chr=="chrM" or len(chr)>5): continue
        begin=int(begin); end=int(end)
        interval=Interval(begin,end)
        interval.chr=chr
        interval.begin=interval.intCenter()-int(CHUNK_WIDTH/2)
        interval.end=interval.intCenter()+int(CHUNK_WIDTH/2)
        interval.type="enhancer"
        interval.id=id
        interval.numPeaks=0
        interval.mates=[]
        interval.peaks=None
        if(hash.get(chr,None) is None): hash[chr]=[]
        hash[chr].append(interval)
    deduplicateEnhancers(hash)
    for chr in hash.keys():
        array=hash[chr]
        for elem in array:
            byChr[chr].append(elem)
            byID[elem.id]=elem
    IN.close()

def deduplicateEnhancers(byChr):
    for chr in byChr.keys():
        array=byChr[chr]
        array.sort(key=lambda x: x.begin)
        L=len(array)
        i=0
        while(i+1<L):
            if(array[i].overlaps(array[i+1])):
                del array[i]
                L-=1
            else: i+=1

def loadStandardized(filename,whichTime,byID):
    with open(filename) as IN:
        for line in IN:
            fields=line.rstrip().split()
            if(len(fields)!=11): continue
            (fastb,LLR,P,parse,features,GR,AP1,CEBP,FOX,KLF,CTCF)=fields
            if(not rex.find("(\S+)\.(t\d+)\.fastb",fastb)):
                raise Exception(fastb)
            enhancerID=rex[1]; time=rex[2]
            if(time!=whichTime): continue
            enhancer=byID.get(enhancerID,None)
            if(enhancer is None): continue
            enhancer.LLR=LLR; enhancer.P=P
            origin=enhancer.begin
            peaks=[]
            fields=parse.split("|")
            for field in fields:
                rex.findOrDie("(\d+):(\d+)-(\d+)",field)
                state=int(rex[1]); begin=int(rex[2]); end=int(rex[3])
                if(state==3): peaks.append(Interval(origin+begin,origin+end))
            fields=features.split("|")
            if(len(fields)!=len(peaks)): raise Exception("unequal")
            for i in range(len(fields)):
                field=fields[i]
                subfields=field.split(",")
                (DNase_t0,DNase_t3,P300_t0,P300_t3)=subfields
                peak=peaks[i]
                peak.std_dnase_t0=round(float(DNase_t0),3)
                peak.std_dnase_t3=round(float(DNase_t3),3)
                peak.std_p300_t0=round(float(P300_t0),3)
                peak.std_p300_t3=round(float(P300_t3),3)
                peak.motifs=set()
            enhancer.peaks=peaks
            enhancer.numPeaks=len(peaks)

def loadRaw(filename,whichTime,byID):
    with open(filename) as IN:
        for line in IN:
            fields=line.rstrip().split()
            if(len(fields)!=11): continue
            (fastb,LLR,P,parse,features,GR,AP1,CEBP,FOX,KLF,CTCF)=fields
            if(not rex.find("(\S+)\.(t\d+)\.fastb",fastb)):
                raise Exception(fastb)
            enhancerID=rex[1]; time=rex[2]
            if(time!=whichTime): continue
            enhancer=byID.get(enhancerID,None)
            if(enhancer is None): continue
            peaks=enhancer.peaks
            fields=features.split("|")
            if(len(fields)!=len(peaks)): raise Exception("unequal")
            for i in range(len(fields)):
                field=fields[i]
                subfields=field.split(",")
                (DNase_t0,DNase_t3,P300_t0,P300_t3)=subfields
                peak=peaks[i]
                peak.raw_dnase_t0=round(float(DNase_t0),3)
                peak.raw_dnase_t3=round(float(DNase_t3),3)
                peak.raw_p300_t0=round(float(P300_t0),3)
                peak.raw_p300_t3=round(float(P300_t3),3)

def sortArrays(byChr):
    chroms=byChr.keys()
    for chr in chroms:
        byChr[chr].sort(key=lambda x: x.begin)

def assignToAnchors(fromType,byChr):
    chroms=byChr.keys()
    for chr in chroms:
        array=byChr[chr]
        L=len(array)
        for i in range(L):
            this=array[i]
            if(this.type!=fromType): continue
            left=getLeftNearest(array,i,"anchor")
            right=getRightNearest(array,i,"anchor")
            if(left and this.begin-left.end<MAX_DIST_TO_ANCHOR):
                left.objects.add(this)
            if(right and right.begin-this.end<MAX_DIST_TO_ANCHOR):
                right.objects.add(this)

def getLeftNearest(array,i,toType):
    origin=array[i].begin
    for j in range(i-1,-1,-1):
        if(array[j].type==toType):
            return array[j]
    return None
    
def getRightNearest(array,i,toType):
    origin=array[i].end
    for j in range(i+1,len(array)):
        if(array[j].type==toType):
            return array[j]
    return None

def buildGraph(byChr):
    sortArrays(byChr)
    assignToAnchors("enhancer",byChr)
    assignToAnchors("gene",byChr)
    pairViaAnchors(byChr)

def pairViaAnchors(byChr):
    for chr in byChr.keys():
        array=byChr[chr]
        for anchor in array:
            if(anchor.type!="anchor"): continue
            for object in anchor.objects:
                if(object.type!="enhancer" and object.type!="gene"): continue
                for other in anchor.mate.objects:
                    object.mates.append(other)

def countAnchorOverlaps(byChr):
    seenGenes=set(); seenEnhancers=set(); seenPairs=set()
    allGenes=set(); allEnhancers=set()
    totalAnchors=0
    for chr in byChr.keys():
        array=byChr[chr]
        for elem in array:
            if(elem.type=="anchor"):
                totalAnchors+=1
                for object in elem.objects:
                    if(object.type=="gene"):
                        seenGenes.add(object.id)
                        for mate in object.mates:
                            if(mate.type=="enhancer"):
                                seenPairs.add(object.id+"_"+mate.id)
                    elif(object.type=="enhancer"):
                        seenEnhancers.add(object.id)
            elif(elem.type=="gene"): allGenes.add(elem.id)
            elif(elem.type=="enhancer"): allEnhancers.add(elem.id)
    overlappedGenes=len(seenGenes)
    overlappedEnhancers=len(seenEnhancers)
    totalGenes=len(allGenes)
    totalEnhancers=len(allEnhancers)
    totalPairs=len(seenPairs)
    print("overlapped enhancers:",overlappedEnhancers,
          "total enhancers:",totalEnhancers,
          "overlapped genes:",overlappedGenes,
          "total genes:",totalGenes,
          "linked pairs:",totalPairs,
          "total anchors:",totalAnchors,
          sep="\t")
    

def loadMotifs(filename,byID):
    IN=open(filename)
    for line in IN:
        fields=line.rstrip().split()
        if(len(fields)<2): continue
        peakID=fields[0]
        motifs=fields[1:]
        rex.findOrDie("(\S+)_(\d+)",peakID)
        enhancerID=rex[1]; peakNum=int(rex[2])
        enhancer=byID.get(enhancerID,None)
        if(enhancer is None): continue
        peak=enhancer.peaks[peakNum]
        for motif in motifs: peak.motifs.add(motif)
    IN.close()

def dumpGraph(byChr):
    chroms=byChr.keys()
    for chr in chroms:
        array=byChr[chr]
        for elem in array:
            if(elem.type=="gene"): dumpGene(elem)
            elif(elem.type=="enhancer"): dumpEnhancer(elem)

def dumpGene(gene):
    essex=makeGeneEssex(gene)
    essex.print(sys.stdout)
    print()

def dumpEnhancer(enhancer):
    essex=makeEnhancerEssex(enhancer)
    if(essex is None): return
    essex.print(sys.stdout)
    print()

def makeGeneEssex(gene):
    essex=EssexNode(["gene"])
    essex.setAttribute("id",gene.id)
    essex.setAttribute("logFC",gene.logFC)
    essex.setAttribute("FDR",gene.FDR)
    essex.setAttribute("chr",gene.chr)
    essex.setAttribute("TSS",gene.begin)
    mates=EssexNode(["mates"])
    essex.addElem(mates)
    for mate in gene.mates:
        mates.addElem(mate.id)
    return essex

def makeEnhancerEssex(enhancer):
    essex=EssexNode(["enhancer"])
    essex.setAttribute("id",enhancer.id)
    essex.setAttribute("chr",enhancer.chr)
    essex.setAttribute("begin",enhancer.begin)
    essex.setAttribute("end",enhancer.end)
    essex.setAttribute("numPeaks",enhancer.numPeaks)
    peaks=EssexNode(["peaks"])
    essex.addElem(peaks)
    if(enhancer.peaks is None): return None
    for peak in enhancer.peaks:
        node=EssexNode(["peak"])
        peaks.addElem(node)
        node.setAttribute("begin",peak.begin)
        node.setAttribute("end",peak.end)
        node.setAttribute("width",peak.end-peak.begin)
        node.setAttribute("std_dnase_t0",peak.std_dnase_t0)
        node.setAttribute("std_dnase_t3",peak.std_dnase_t3)
        node.setAttribute("raw_dnase_t0",peak.raw_dnase_t0)
        node.setAttribute("raw_dnase_t3",peak.raw_dnase_t3)
        node.setAttribute("std_P300_t0",peak.std_p300_t0)
        node.setAttribute("std_P300_t3",peak.std_p300_t3)
        node.setAttribute("raw_P300_t0",peak.raw_p300_t0)
        node.setAttribute("raw_P300_t3",peak.raw_p300_t3)
        motifsNode=EssexNode(["motifs"])
        node.addElem(motifsNode)
        for motif in peak.motifs: motifsNode.addElem(motif)
    mates=EssexNode(["mates"])
    essex.addElem(mates)
    for mate in enhancer.mates:
        mates.addElem(mate.id)
    return essex

#=========================================================================
# main()
#=========================================================================
if(len(sys.argv)!=9):
    exit(ProgramName.get()+" <hi-c.bedpe> <tss.txt> <edgeR.txt> <enhancers.bed> <genomewide-features-standardized.txt> <genomewide-features-raw.bed> <motifs-in-peaks.txt> <timepoint>\n")
(hiCfile,tssFile,edgeRfile,enhancerFile,featuresFileStandard,featuresFileRaw,
 motifFile,timepoint)=sys.argv[1:]

byChr=loadLoops(hiCfile)
expressionHash=loadExpressionFile(edgeRfile)
loadTssFile(tssFile,byChr,expressionHash)
byEnhancerID={}
loadEnhancers(enhancerFile,byChr,byEnhancerID)
#loadStandardized(featuresFileStandard,timepoint,byEnhancerID)
#loadRaw(featuresFileRaw,timepoint,byEnhancerID)
#loadMotifs(motifFile,byEnhancerID)
buildGraph(byChr)
countAnchorOverlaps(byChr)
#dumpGraph(byChr)



