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

MAX_DIST_TO_ANCHOR=1
MAX_LINEAR_DIST=10000
#MAX_LINEAR_DIST=35000
MISSING=1000000000
SEEN=set()

def loadLoops(filename):
    byChr={}
    IN=open(filename,"rt")
    for line in IN:
        fields=line.rstrip().split()
        if(len(fields)!=6): continue
        (fromChr,fromBegin,fromEnd,toChr,toBegin,toEnd)=fields
        fromBegin=int(fromBegin); fromEnd=int(fromEnd)
        toBegin=int(toBegin); toEnd=int(toEnd)
        interval1=makeLoopAnchor(fromChr,fromBegin,fromEnd,toChr,toBegin,toEnd)
        interval2=makeLoopAnchor(toChr,toBegin,toEnd,fromChr,fromBegin,fromEnd)
        interval1.mate=interval2
        interval2.mate=interval1
        addAnchor(interval1,byChr)
        addAnchor(interval2,byChr)
    IN.close()
    return byChr

def addAnchor(interval,hash):
    if(hash.get(interval.chr,None) is None): hash[interval.chr]=[]
    hash[interval.chr].append(interval)

def makeLoopAnchor(fromChr,fromBegin,fromEnd,toChr,toBegin,toEnd):
    interval=Interval(fromBegin,fromEnd)
    interval.chr=fromChr
    interval.otherChr=toChr
    interval.otherBegin=toBegin
    interval.otherEnd=toEnd
    interval.type="anchor"
    interval.objects=set()
    return interval

def loadTssFile(filename,byChr):
    IN=open(filename,"rt")
    for line in IN:
        fields=line.rstrip().split()
        if(len(fields)!=4): continue
        (chr,pos,strand,id)=fields
        if(chr=="chrM"): continue
        pos=int(pos)
        interval=Interval(pos,pos+1)
        interval.type="gene"
        interval.id=id
        byChr[chr].append(interval)        
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

def loadEnhancers(filename,byChr):
    IN=open(filename,"rt")
    hash={}
    for line in IN:
        fields=line.rstrip().split()
        if(len(fields)!=4): continue
        (chr,begin,end,id)=fields
        if(chr=="chrM" or len(chr)>5): continue
        begin=int(begin); end=int(end)
        interval=Interval(begin,end)
        interval.type="enhancer"
        interval.id=id
        if(hash.get(chr,None) is None): hash[chr]=[]
        hash[chr].append(interval)
    deduplicateEnhancers(hash)
    for chr in hash.keys():
        array=hash[chr]
        for elem in array:
            byChr[chr].append(elem)
    IN.close()

def distanceToNearest(fromType,toType,byChr):
    chroms=byChr.keys()
    distances=[]
    for chr in chroms:
        array=byChr[chr]
        L=len(array)
        for i in range(L):
            this=array[i]
            if(this.type!=fromType): continue
            leftDist=getLeftDistance(array,i,toType)
            rightDist=getRightDistance(array,i,toType)
            distance=leftDist if leftDist<rightDist else rightDist
            if(distance==MISSING): continue
            distances.append(distance)
    writeValues(distances,"distances.tmp")

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

def writeValues(array,filename):
    OUT=open(filename,"wt")
    for x in array: print(x,file=OUT)
    OUT.close()

def getLeftDistance(array,i,toType):
    origin=array[i].begin
    for j in range(i-1,-1,-1):
        if(array[j].type==toType):
            d=origin-array[j].end
            return d if d>0 else 0
    return MISSING
    
def getRightDistance(array,i,toType):
    origin=array[i].end
    for j in range(i+1,len(array)):
        if(array[j].type==toType):
            d=array[j].begin-origin
            return d if d>0 else 0
    return MISSING
    
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
    
def sortArrays(byChr):
    chroms=byChr.keys()
    for chr in chroms:
        byChr[chr].sort(key=lambda x: x.begin)

def emitPairs(byChr):
    chroms=byChr.keys()
    for chr in chroms:
        array=byChr[chr]
        for elem in array:
            if(elem.type=="anchor"): emitAnchorPairs(elem)

def getType(objects,type):
    array=[]
    for elem in objects:
        if(elem.type==type): array.append(elem)
    return array

def emitGenesEnhancers(objects1,objects2):
    genes=getType(objects1,"gene")
    enhancers=getType(objects2,"enhancer")
    for gene in genes:
        for enhancer in enhancers:
            print(gene.id,enhancer.id,sep="\t")

def emitAnchorPairs(anchor):
    emitGenesEnhancers(anchor.objects,anchor.mate.objects)
    #emitGenesEnhancers(anchor.mate.objects,anchor.objects)

def emitByProximity(byChr):
    global SEEN
    for chr in byChr.keys():
        array=byChr[chr]
        for i in range(len(array)):
            if(array[i].type!="gene"): continue
            gene=array[i]
            if(gene.id in SEEN): continue
            left=getLeftNearest(array,i,"enhancer")
            right=getRightNearest(array,i,"enhancer")
            nearest=left
            if(not left or right and 
               right.distance(gene)<left.distance(gene)):
                nearest=right
            if(nearest and nearest.distance(gene)<MAX_LINEAR_DIST):
                print(gene.id,nearest.id,sep="\t")
                SEEN.add(gene.id)

def countPairs(type1,type2,byChr):
    pairs=set()
    for chr in byChr.keys():
        array=byChr[chr]
        for anchor in array:
            if(anchor.type!="anchor"): continue
            getPairKeys(anchor.objects,anchor.mate.objects,type1,type2,pairs)
            #getPairKeys(anchor.mate.objects,anchor.objects,type1,type2,pairs)
    return len(pairs)

def getPairKeys(objects1,objects2,type1,type2,pairs):
    objects1=getType(objects1,type1)
    objects2=getType(objects2,type2)
    for a in objects1:
        for b in objects2:
            if(a==b): continue
            if(type1==type2):
                if(a.id>b.id): c=a; a=b; b=c
            key=a.id+" "+b.id
            pairs.add(key)
                
#=========================================================================
# main()
#=========================================================================
if(len(sys.argv)!=4):
    exit(ProgramName.get()+" <enhancers.bed> <tss.txt> <loops.txt>\n")
(enhancerFile,tssFile,loopFile)=sys.argv[1:]

byChr=loadLoops(loopFile)
loadTssFile(tssFile,byChr)
loadEnhancers(enhancerFile,byChr)
sortArrays(byChr)
assignToAnchors("enhancer",byChr)
assignToAnchors("gene",byChr)
#emitPairs(byChr)
#emitByProximity(byChr)

#print("gene-gene:",countPairs("gene","gene",byChr))
#print("enhancer-enhancer:",countPairs("enhancer","enhancer",byChr))
#print("enhancer-gene:",countPairs("enhancer","gene",byChr))
#print("gene-enhancer:",countPairs("gene","enhancer",byChr))


