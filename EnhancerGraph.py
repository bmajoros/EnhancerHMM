#=========================================================================
# This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
# License (GPL) version 3, as described at www.opensource.org.
# Copyright (C)2016 William H. Majoros (martiandna@gmail.com).
#=========================================================================
from __future__ import (absolute_import, division, print_function,
   unicode_literals, generators, nested_scopes, with_statement)
from builtins import (bytes, dict, int, list, object, range, str, ascii,
   chr, hex, input, next, oct, open, pow, round, super, filter, map, zip)
import sys
import ProgramName
from Interval import Interval
from EssexNode import EssexNode
from Rex import Rex

#=========================================================================
# Attributes:
#   byChr = hash : chromosome to Intervals with attribute type indicating
#                  "enhancer", "gene", or "anchor"
# Instance Methods:
#   graph=EnhancerGraph(enhancerTimepoint,geneTimepoint)
#   genes=graph.getGenes()
#   enhancers=graph.getEnhancers()
#   singletons=graph.getSingletons()
#   multipeakEnhancers=graph.getMultipeakEnhancers()
#   byChr=graph.getByChr() # hash : chromosome => array of genes, enhancers,
#                                                 and loop anchors
#   hash=graph.getEnhancerHash() # hash : enhancerID => array of enhancers
#   TFs=graph.getAllTFs()
#
# Class Methods:
#   
#=========================================================================
class EnhancerGraph:
    """EnhancerGraph"""
    def __init__(self,enhancerTimepoint,geneTimepoint):
        self.rex=Rex()
        self.CHUNK_WIDTH=1000
        self.MAX_DIST_TO_ANCHOR=1
        #hiCfile="/data/reddylab/Tony/HiC/GC_Timecourse/merged_hics/iter3/merged.loops.fdr.0.01.bedpe"
        hiCfile="/data/reddylab/Tony/HiC/GC_Timecourse/merged_hics/iter3/merged.loops.fdr.0.001.bedpe"
        tssFile="tss.hg38.txt"
        edgeRfile="/data/reddylab/projects/GGR/results/rna_seq/differential_expression/iter0/edgeR/edgeR.sva."+geneTimepoint+".vs.t00.protein_coding.txt"
        enhancerFile="p300-hg38.bed"
        featuresFileStandard="genomewide-features-standardized.txt"
        featuresFileRaw="genomewide-features-raw.txt"
        motifFile="motifs-in-peaks-"+enhancerTimepoint+".txt"
        self.byChr=self.loadLoops(hiCfile)
        expressionHash=self.loadExpressionFile(edgeRfile)
        self.loadTssFile(tssFile,self.byChr,expressionHash)
        self.byEnhancerID={}
        self.loadEnhancers(enhancerFile,self.byChr,self.byEnhancerID)
        self.loadStandardized(featuresFileStandard,enhancerTimepoint,
                              self.byEnhancerID)
        self.loadRaw(featuresFileRaw,enhancerTimepoint,self.byEnhancerID)
        self.loadMotifs(motifFile,self.byEnhancerID)
        self.getGenesAndEnhancers()
        self.buildGraph(self.byChr)
        #dumpGraph(byChr)

    def getSingletons(self):
        array=[]
        for enhancer in self.enhancers:
            if(enhancer.numPeaks==1): array.append(enhancer)
        return array

    def getMultipeakEnhancers(self):
        array=[]
        for enhancer in self.enhancers:
            if(enhancer.numPeaks>1): array.append(enhancer)
        return array

    def getAllTFs(self):
        all=set()
        for enhancer in self.enhancers:
            for peak in enhancer.peaks:
                for motif in peak.motifs:
                    all.add(motif)
        return all

    def getGenes(self):
        return self.genes

    def getEnhancers(self):
        return self.enhancers

    def getGenesAndEnhancers(self):
        self.genes=[]; self.enhancers=[]
        for chr in self.byChr.keys():
            array=self.byChr[chr]
            for elem in array:
                if(elem.type=="gene"): self.genes.append(elem)
                elif(elem.type=="enhancer"): self.enhancers.append(elem)

    def getByChr(self):
        return self.byChr

    def getEnhancerHash(self):
        return self.byEnhancerID
    
    def loadLoops(self,filename):
        byChr={}
        IN=open(filename,"rt")
        for line in IN:
            fields=line.rstrip().split()
            if(len(fields)!=6): continue
            (fromChr,fromBegin,fromEnd,toChr,toBegin,toEnd)=fields
            fromBegin=int(fromBegin); fromEnd=int(fromEnd)
            toBegin=int(toBegin); toEnd=int(toEnd)
            interval1=self.makeLoopAnchor(fromChr,fromBegin,fromEnd)
            interval2=self.makeLoopAnchor(toChr,toBegin,toEnd)
            interval1.mate=interval2
            interval2.mate=interval1
            self.addAnchor(interval1,byChr)
            self.addAnchor(interval2,byChr)
        IN.close()
        return byChr
    
    def addAnchor(self,interval,hash):
        if(hash.get(interval.chr,None) is None): hash[interval.chr]=[]
        hash[interval.chr].append(interval)
    
    def makeLoopAnchor(self,fromChr,fromBegin,fromEnd):
        interval=Interval(fromBegin,fromEnd)
        interval.chr=fromChr
        interval.type="anchor"
        interval.objects=set()
        return interval
    
    def loadExpressionFile(self,filename):
        hash={}
        with open(filename) as IN:
            for line in IN:
                fields=line.rstrip().split()
                if(len(fields)!=6): continue
                (gene,logFC,logCPM,LR,Pvalue,FDR)=fields
                hash[gene]=[logFC,logCPM,LR,Pvalue,FDR]
        return hash
    
    def loadTssFile(self,filename,byChr,expressionHash):
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
    
    def loadEnhancers(self,filename,byChr,byID):
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
            interval.begin=interval.intCenter()-int(self.CHUNK_WIDTH/2)
            interval.end=interval.intCenter()+int(self.CHUNK_WIDTH/2)
            interval.type="enhancer"
            interval.id=id
            interval.numPeaks=0
            interval.mates=[]
            interval.peaks=[]
            if(hash.get(chr,None) is None): hash[chr]=[]
            hash[chr].append(interval)
        self.deduplicateEnhancers(hash)
        for chr in hash.keys():
            array=hash[chr]
            for elem in array:
                byChr[chr].append(elem)
                byID[elem.id]=elem
        IN.close()
    
    def deduplicateEnhancers(self,byChr):
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
    
    def loadStandardized(self,filename,whichTime,byID):
        rex=self.rex
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
                    if(state==3):
                        peaks.append(Interval(origin+begin,origin+end))
                fields=features.split("|")
                if(len(fields)!=len(peaks)): raise Exception("unequal")
                for i in range(len(fields)):
                    field=fields[i]
                    subfields=field.split(",")
                    (DNase_t0,DNase_t3,P300_t0,P300_t3)=subfields
                    peak=peaks[i]
                    peak.std_dnase_t0=float(DNase_t0)
                    peak.std_dnase_t3=float(DNase_t3)
                    peak.std_p300_t0=float(P300_t0)
                    peak.std_p300_t3=float(P300_t3)
                    peak.motifs=set()
                    peak.open=peak.std_dnase_t3>=1.0
                    peak.active=peak.std_p300_t3>=1.0
                    peak.dex_responsive=peak.active and peak.std_p300_t0<1.0
                enhancer.peaks=peaks
                enhancer.numPeaks=len(peaks)
    
    def loadRaw(self,filename,whichTime,byID):
        rex=self.rex
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
                    peak.raw_dnase_t0=float(DNase_t0)
                    peak.raw_dnase_t3=float(DNase_t3)
                    peak.raw_p300_t0=float(P300_t0)
                    peak.raw_p300_t3=float(P300_t3)
    
    def sortArrays(self,byChr):
        chroms=byChr.keys()
        for chr in chroms:
            byChr[chr].sort(key=lambda x: x.begin)
    
    def assignToAnchors(self,fromType,byChr):
        chroms=byChr.keys()
        for chr in chroms:
            array=byChr[chr]
            L=len(array)
            for i in range(L):
                this=array[i]
                if(this.type!=fromType): continue
                left=self.getLeftNearest(array,i,"anchor")
                right=self.getRightNearest(array,i,"anchor")
                if(left and this.begin-left.end<self.MAX_DIST_TO_ANCHOR):
                    left.objects.add(this)
                if(right and right.begin-this.end<self.MAX_DIST_TO_ANCHOR):
                    right.objects.add(this)
    
    def getLeftNearest(self,array,i,toType):
        origin=array[i].begin
        for j in range(i-1,-1,-1):
            if(array[j].type==toType):
                return array[j]
        return None
        
    def getRightNearest(self,array,i,toType):
        origin=array[i].end
        for j in range(i+1,len(array)):
            if(array[j].type==toType):
                return array[j]
        return None
    
    def buildGraph(self,byChr):
        self.sortArrays(byChr)
        self.assignToAnchors("enhancer",byChr)
        self.assignToAnchors("gene",byChr)
        self.pairViaAnchors(byChr)
        self.classifyMateTypes(byChr)

    def classifyMateTypes(self,byChr):
        isEnhancer=set(); isGene=set()
        for enhancer in self.enhancers: isEnhancer.add(enhancer.id)
        for gene in self.genes: isGene.add(gene.id)
        for chr in byChr.keys():
            array=byChr[chr]
            for elem in array:
                if(elem.type=="anchor"): continue
                elem.geneMates=[]
                elem.enhancerMates=[]
                for mate in elem.mates:
                    if(mate.id in isEnhancer): elem.enhancerMates.append(mate)
                    elif(mate.id in isGene): elem.geneMates.append(mate)
    
    def pairViaAnchors(self,byChr):
        for chr in byChr.keys():
            array=byChr[chr]
            for anchor in array:
                if(anchor.type!="anchor"): continue
                for object in anchor.objects:
                    if(object.type!="enhancer" and object.type!="gene"):
                        continue
                    for other in anchor.mate.objects:
                        object.mates.append(other)
    
    def loadMotifs(self,filename,byID):
        rex=self.rex
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
            for motif in motifs:
                m=motif
                if(m=="NR3C1"): m="GR"
                peak.motifs.add(m)
        IN.close()
    
    def dumpGraph(self,byChr):
        chroms=byChr.keys()
        for chr in chroms:
            array=byChr[chr]
            for elem in array:
                if(elem.type=="gene"): self.dumpGene(elem)
                elif(elem.type=="enhancer"): self.dumpEnhancer(elem)
    
    def dumpGene(self,gene):
        essex=self.makeGeneEssex(gene)
        essex.print(sys.stdout)
        print()
    
    def dumpEnhancer(self,enhancer):
        essex=self.makeEnhancerEssex(enhancer)
        if(essex is None): return
        essex.print(sys.stdout)
        print()
    
    def makeGeneEssex(self,gene):
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
    
    def makeEnhancerEssex(self,enhancer):
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
    


