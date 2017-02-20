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
from scipy.stats import mannwhitneyu
from Interval import Interval
from SummaryStats import SummaryStats
from Rex import Rex
rex=Rex()

P300_TYPE="change" # "change" or "first-time" or "second-time"
P300hash={}
MIN_PEAK_LEN=1 # 250
MAX_DIST_TO_ANCHOR=1
MAX_LINEAR_DIST=10000
#MAX_LINEAR_DIST=35000
MISSING=1000000000
SEEN=set()

def getNumPeaks(filename,whichTime):
    global P300hash
    # iter0_peak14413.t3.fastb        2234.2906200000034      1.0     1:0-1|2:1-2|3:2-1243|4:1243-1256|3:1256-1998|4:1998-1999|5:1999-2000    1.07942122577,1.12211048111,1.00477844983,1.16141454506|0.952633976593,1.02660980228,0.892990209161,0.961271510352      10      10      11      10      10      00
    hash={} # maps enhancerID to P300 value
    with open(filename) as IN:
        for line in IN:
            fields=line.rstrip().split()
            if(len(fields)!=11): continue
            (fastb,LLR,P,parse,features,GR,AP1,CEBP,FOX,KLF,CTCF)=fields
            ### if(float(LLR)<1000): continue
            if(not rex.find("(\S+)\.(t\d+)\.fastb",fastb)):
                raise Exception(fastb)
            enhancerID=rex[1]
            time=rex[2]
            if(time!=whichTime): continue
            peakLengths=[]
            fields=parse.split("|")
            for field in fields:
                if(not rex.find("(\d+):(\d+)-(\d+)",field)): exit(field)
                state=int(rex[1]); begin=int(rex[2]); end=int(rex[3])
                if(state==3): peakLengths.append(end-begin)
            #numPeaks=len(fields)
            numPeaks=0
            for L in peakLengths:
                if(L>=MIN_PEAK_LEN): numPeaks+=1
            hash[enhancerID]=numPeaks

            fields=features.split("|")
            array=[]
            for i in range(len(fields)):
                if(peakLengths[i]<MIN_PEAK_LEN): continue
                field=fields[i]
                subfields=field.split(",")
                (DNase_t0,DNase_t3,P300_t0,P300_t3)=subfields
                delta=float(P300_t3)-float(P300_t0)
                if(P300_TYPE=="change"): array.append(delta)
                elif(P300_TYPE=="first-time"): array.append(float(P300_t0))
                elif(P300_TYPE=="second-time"): array.append(float(P300_t3))
                else: raise Exception(P300_TYPE)
            if(len(array)>0):
                P300hash[enhancerID]=max(array)
                #P300hash[enhancerID]=sum(array)
    return hash

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
        interval.linkedTo=[]
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
        interval.begin=interval.intCenter()-1000
        interval.end=interval.intCenter()+1000
        interval.type="enhancer"
        interval.id=id
        interval.numPeaks=-1
        interval.linkedTo=[]
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

def pairWithEnhancers(byChr,numPeaksHash):
    for chr in byChr.keys():
        array=byChr[chr]
        for anchor in array:
            if(anchor.type!="anchor"): continue
            for enhancer in anchor.objects:
                if(enhancer.type!="enhancer"): continue
                if(numPeaksHash.get(enhancer.id,None) is None):
                    enhancer.numPeaks=0
                else: enhancer.numPeaks=numPeaksHash[enhancer.id]
                for object in anchor.mate.objects:
                    enhancer.linkedTo.append(object)
            for gene in anchor.objects:
                if(gene.type!="gene"): continue
                for object in anchor.mate.objects:
                    gene.linkedTo.append(object)

def getDumpType(object):
    if(object.type=="enhancer"):
        if(object.numPeaks<1): return "bad"
        elif(object.numPeaks==1): return "singleton"
        return "multipeak"
    return object.type

def dumpPairedP300(byChr):
    MULTIMULTI=open("multimulti.txt","wt")
    MULTISINGLE=open("multisingle.txt","wt")
    SINGLESINGLE=open("singlesingle.txt","wt")
    for chr in byChr.keys():
        for enhancer in byChr[chr]:
            if(enhancer.type!="enhancer" or enhancer.numPeaks<1): continue
            if(P300hash.get(enhancer.id,None) is None): continue
            p300_1=P300hash[enhancer.id]
            fromType=getDumpType(enhancer)
            mates=enhancer.linkedTo
            for mate in mates:
                if(mate.type!="enhancer" or mate.numPeaks<1): continue
                toType=getDumpType(mate)
                if(fromType==toType and mate.begin<enhancer.begin): continue
                if(toType=="bad"): continue
                if(P300hash.get(mate.id,None) is None): continue
                p300_2=P300hash[mate.id]
                if(fromType=="multipeak" and toType=="multipeak"):
                    print(p300_1,p300_2,sep="\t",file=MULTIMULTI)
                if(fromType=="multipeak" and toType=="singleton"):
                    print(p300_1,p300_2,sep="\t",file=MULTISINGLE)
                if(fromType=="singleton" and toType=="singleton"):
                    print(p300_1,p300_2,sep="\t",file=SINGLESINGLE)
    MULTIMULTI.close(); MULTISINGLE.close(); SINGLESINGLE.close()

def dumpLinks(byChr):
    numSingletons=0
    numMultipeaks=0
    numGenes=0
    counts={}
    degreeByType={}
    clusterSeen=set()
    clusterCounts={}
    CONNECTED=open("connected-to-multi.txt","wt")
    NOT_CONNECTED=open("not-connected-to-multi.txt","wt")
    for chr in byChr.keys():
        for object in byChr[chr]: getClusters(object,clusterSeen,clusterCounts)
        for enhancer in byChr[chr]:
            if(enhancer.type!="enhancer" or enhancer.numPeaks<1): continue
            if(enhancer.numPeaks==1): numSingletons+=1
            if(enhancer.numPeaks>1): numMultipeaks+=1
            fromType=getDumpType(enhancer)
            mates=enhancer.linkedTo
            degree=0
            connectedToMulti=False
            for mate in mates:
                toType=getDumpType(mate)
                if(toType=="bad"): continue
                if(fromType==toType and mate.begin<enhancer.begin): continue
                if(counts.get(fromType,None) is None):
                    counts[fromType]={}
                counts[fromType][toType]=counts[fromType].get(toType,0)+1
                degree+=1
                if(toType=="multipeak"): connectedToMulti=True
            if(fromType=="singleton"):
                FILE=CONNECTED if connectedToMulti else NOT_CONNECTED
                print(enhancer.id,file=FILE)
            if(degreeByType.get(fromType,None) is None):
                degreeByType[fromType]=[]
            degreeByType[fromType].append(degree)
        for gene in byChr[chr]:
            if(gene.type!="gene"): continue
            fromType="gene"
            numGenes+=1
            mates=gene.linkedTo
            degree=0
            for mate in mates:
                toType=getDumpType(mate)
                if(toType=="bad"): continue
                if(counts.get(fromType,None) is None):
                    counts[fromType]={}
                counts[fromType][toType]=counts[fromType].get(toType,0)+1
                degree+=1
            if(degreeByType.get(fromType,None) is None):
                degreeByType[fromType]=[]
            degreeByType[fromType].append(degree)
    print(numSingletons,"total singletons")
    print(numMultipeaks,"total multipeaks")
    print(numGenes,"total genes")
    fromTypes=counts.keys()
    for fromType in fromTypes:
        (mean,SD,min,max)=\
            SummaryStats.roundedSummaryStats(degreeByType[fromType])
        print("degree",fromType,mean,"+/-",SD)
    mannWhitney(degreeByType["singleton"],degreeByType["multipeak"],
                "MannWhitney: singleton-multipeak")
    mannWhitney(degreeByType["gene"],degreeByType["multipeak"],
                "MannWhitney: gene-multipeak")
    mannWhitney(degreeByType["gene"],degreeByType["singleton"],
                "MannWhitney: gene-singleton")
    for fromType in fromTypes:
        hash=counts[fromType]
        toTypes=hash.keys()
        for toType in toTypes:
            count=hash[toType]
            print(fromType,toType,count,sep="\t")
    for key in clusterCounts.keys():
        value=clusterCounts[key]
        print("cluster",key,value)

def mannWhitney(array1,array2,label):
    (U,P)=mannwhitneyu(array1,array2)
    print(label,"P="+str(P),"U="+str(U),sep="\t")

def getClusters(object,seen,counts):
    if(object.type=="anchor"): return
    if(object.id in seen): return
    cluster=set()
    stack=[]
    stack.append(object)
    while(len(stack)>0):
        last=len(stack)-1
        object=stack[last]
        del stack[last]
        seen.add(object.id)
        t=getDumpType(object)
        if(t=="bad"): continue
        cluster.add(t)
        for mate in object.linkedTo:
            if(mate.id in seen): continue
            if(getDumpType(mate)=="bad"): continue
            stack.append(mate)
    array=[]
    for elem in cluster: array.append(elem)
    array.sort()
    key=""
    for elem in array: key+=elem+" "
    counts[key]=counts.get(key,0)+1

#=========================================================================
# main()
#=========================================================================
if(len(sys.argv)!=6):
    exit(ProgramName.get()+" <genomewide-features.txt> <enhancers.bed> <tss.txt> <loops.txt> <t#>\n")
(featuresFile,enhancerFile,tssFile,loopFile,timepoint)=sys.argv[1:]

#numPeaksHash=getNumPeaks(featuresFile,"t05")

numPeaksHash=getNumPeaks(featuresFile,timepoint) ###

byChr=loadLoops(loopFile)
loadTssFile(tssFile,byChr)
loadEnhancers(enhancerFile,byChr)
sortArrays(byChr)
assignToAnchors("enhancer",byChr)
assignToAnchors("gene",byChr)
pairWithEnhancers(byChr,numPeaksHash)
dumpLinks(byChr)
dumpPairedP300(byChr)

#print("gene-gene:",countPairs("gene","gene",byChr))
#print("enhancer-enhancer:",countPairs("enhancer","enhancer",byChr))
#print("enhancer-gene:",countPairs("enhancer","gene",byChr))
#print("gene-enhancer:",countPairs("gene","enhancer",byChr))


