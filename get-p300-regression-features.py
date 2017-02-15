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
from Shuffler import Shuffler
from SummaryStats import SummaryStats
from Pipe import Pipe
from Rex import Rex
rex=Rex()

MIN_PEAK_LEN=250
WANT_SHUFFLE=False
MIN_DNASE=1.0 # standardized
MIN_P300=1.0  # standardized
PEAK_STATE=3


class Enhancer:
    def __init__(self,substrate,peaks):
        self.substrate=substrate
        self.peaks=peaks
    def hasInactive(self):
        for peak in self.peaks:
            if(not peak.active): return True
        return False
    def hasActive(self):
        for peak in self.peaks:
            if(peak.active): return True
        return False
    def getResponsivePeaks(self):
        ret=[]
        for peak in self.peaks:
            if(peak.dex_responsive): ret.append(peak)
        return ret
    def findResponsiveWithMotif(self,motif):
        found=[]
        for peak in self.peaks:
            if(peak.dex_responsive and motif in peak.motifs):
                found.append(peak)
        return found
    def hasPeakWithMotif(self,motif):
        for peak in self.peaks:
            if(motif in peak.motifs): return True
        return False
    def countPeaksWithMotif(self,motif):
        count=0
        for peak in self.peaks:
            if(motif in peak.motifs): count+=1
        return count
    def countOtherPeaksWithMotif(self,motif,ignorePeak):
        count=0
        for peak in self.peaks:
            if(peak is not ignorePeak and
               motif in peak.motifs): count+=1
        return count
    def numActivePeaks(self):
        count=0
        for peak in self.peaks:
            if(peak.active): count+=1
        return count
    def numInactivePeaks(self):
        count=0
        for peak in self.peaks:
            if(not peak.active): count+=1
        return count
    def hasOpen(self):
        for peak in self.peaks:
           if(peak.open): return True
        return False

def load(filename):
    activeButClosed=0
    active=0
    enhancers=[]
    IN=open(filename,"rt")
    for line in IN:
        fields=line.rstrip().split()
        if(len(fields)!=11): continue
        (substrate,LLR,P,parse,features,GR,AP1,CEBP,FOX,KLF,CTCF)=fields
        if(float(P)<0.9): continue
        if(not rex.find("\.t05\.",substrate)): continue
        parseElems=parse.split("|")
        peaks=[]
        for parseElem in parseElems:
            if(not rex.find("(\d+):(\d+)-(\d+)",parseElem)):
                raise Exception(parseElem)
            state=int(rex[1])
            if(state!=PEAK_STATE): continue
            begin=int(rex[2]); end=int(rex[3])
            peaks.append(Interval(begin,end))
        featureSets=features.split("|")
        if(len(featureSets)!=len(peaks)): raise Exception("unequal")
        N=len(featureSets)
        for i in range(N):
            featureSet=featureSets[i]
            feats=featureSet.split(",")
            if(len(feats)!=4): raise Exception(featureSet)
            peak=peaks[i]
            peak.dnase_0=float(feats[0])
            peak.dnase_3=float(feats[1])
            peak.p300_0=float(feats[2])
            peak.p300_3=float(feats[3])
            peak.open=peak.dnase_3>=MIN_DNASE
            peak.active=peak.p300_3>=MIN_P300
            peak.dex_responsive=peak.active and peak.p300_0<MIN_P300
            peak.nearInactive=False
            peak.nearActive=False
            peak.nearResponsive=False
            peak.nearNonResponsive=False
            if(peak.active and not peak.open): activeButClosed+=1
            if(peak.active): active+=1
            peak.motifs=set()
            if(GR[i]=="1"): peak.motifs.add("GR")
            if(AP1[i]=="1"): peak.motifs.add("AP1")
            if(CEBP[i]=="1"): peak.motifs.add("CEBP")
            if(FOX[i]=="1"): peak.motifs.add("FOX")
            if(KLF[i]=="1"): peak.motifs.add("KLF")
            if(CTCF[i]=="1"): peak.motifs.add("CTCF")
        filtered=[]
        for peak in peaks:
            if(peak.length()>=MIN_PEAK_LEN): filtered.append(peak)
        peaks=filtered
        enhancers.append(Enhancer(substrate,peaks))
        N=len(peaks)
        for i in range(N):
            peak=peaks[i]
            if(peak.active):
                if(i>0): peaks[i-1].nearActive=True
                if(i+1<N):peaks[i+1].nearActive=True
            else:
                if(i>0): peaks[i-1].nearInactive=True
                if(i+1<N):peaks[i+1].nearInactive=True
            if(peak.dex_responsive):
                if(i>0): peaks[i-1].nearResponsive=True
                if(i+1<N):peaks[i+1].nearResponsive=True
            else:
                if(i>0): peaks[i-1].nearNonResponsive=True
                if(i+1<N):peaks[i+1].nearNonResponsive=True
    IN.close()
    return enhancers

def getMultipeaks(enhancers):
    multipeaks=[]
    for enhancer in enhancers:
        if(len(enhancer.peaks)>1): multipeaks.append(enhancer)
    return multipeaks

def getSingletons(enhancers):
    singletons=[]
    for enhancer in enhancers:
        if(len(enhancer.peaks)==1): singletons.append(enhancer)
    return singletons

def shuffle(enhancers):
    pool=[]
    for enhancer in enhancers: pool.extend(enhancer.peaks)
    Shuffler.shuffle(pool)
    next=0
    for enhancer in enhancers:
        n=len(enhancer.peaks)
        enhancer.peaks=[]
        for i in range(n):
            enhancer.peaks.append(pool[next])
            next+=1

def getPvalue(array,x):
    count=0
    for value in array:
        if(value>=x): count+=1
    P=float(count)/float(len(array))
    return (P,count)

def mergeData(enhancers,raw):
    results=[]
    hash={}
    for enhancer in raw: hash[enhancer.substrate]=enhancer
    for enhancer in enhancers:
        raw=hash.get(enhancer.substrate,None)
        if(raw is None): continue
        n1=len(enhancer.peaks)
        n2=len(raw.peaks)
        if(n1!=n2): exit(str(n1)+"!="+str(n2))
        for i in range(n1):
            enhancer.peaks[i].rawP300_0=raw.peaks[i].p300_0
            enhancer.peaks[i].rawP300_3=raw.peaks[i].p300_3
            enhancer.peaks[i].deltaP300=\
                enhancer.peaks[i].rawP300_3-enhancer.peaks[i].rawP300_0
            enhancer.peaks[i].rawDNase_0=raw.peaks[i].dnase_0
            enhancer.peaks[i].rawDNase_3=raw.peaks[i].dnase_3
            enhancer.peaks[i].deltaDNase=\
                enhancer.peaks[i].rawDNase_3-enhancer.peaks[i].rawDNase_0
        results.append(enhancer)
    return results

def writeValues(values,filename):
    OUT=open(filename,"wt")
    for value in values:
        OUT.write(str(value)+"\n")
    OUT.close()

def extractFeatures(enhancer):
    numPeaks=len(enhancer.peaks)
    for i in range(numPeaks):
        peak=enhancer.peaks[i]
        #if(not peak.dex_responsive): continue
        features=[]
        features.append(peak.deltaP300)
        features.append(peak.deltaDNase)
        features.append(1 if "GR" in peak.motifs else 0)
        features.append(1 if "AP1" in peak.motifs else 0)
        features.append(1 if "CEBP" in peak.motifs else 0)
        features.append(1 if "FOX" in peak.motifs else 0)
        features.append(1 if "KLF" in peak.motifs else 0)
        features.append(1 if "CTCF" in peak.motifs else 0)
        features.append(enhancer.countOtherPeaksWithMotif("GR",peak))
        features.append(enhancer.countOtherPeaksWithMotif("AP1",peak))
        features.append(enhancer.countOtherPeaksWithMotif("CEBP",peak))
        features.append(enhancer.countOtherPeaksWithMotif("FOX",peak))
        features.append(enhancer.countOtherPeaksWithMotif("KLF",peak))
        features.append(enhancer.countOtherPeaksWithMotif("CTCF",peak))
        nonRespNear=0
        for j in range(numPeaks):
            if(j!=i and not enhancer.peaks[j].dex_responsive): nonRespNear+=1
        features.append(nonRespNear)
        features.append(numPeaks)
        for i in range(len(features)):
            feature=features[i]
            print(feature,end="")
            if(i<len(features)-1): print("\t",end="")
        print()

# dP300 = dDNase
#        + #GRE_here + #GRE_nearby + #AP1_nearby + #NonResponsiveNearby
#        + #CEBP_here + #FOX_here + #KLF_here + #CTCF_here


#=========================================================================
# main()
#=========================================================================
if(len(sys.argv)!=3):
    exit(ProgramName.get()+" <features-standardized.txt> <features-raw.txt>\n")
(standardized,raw)=sys.argv[1:]

# Load data
enhancers=load(standardized)
rawEnhancers=load(raw)
enhancers=mergeData(enhancers,rawEnhancers)

# Extract features
print("dP300\tdDNase\tGR\tAP1\tCEBP\tFOX\tKLF\tCTCF\tGR_nearby\tAP1_nearby\tCEBP_nearby\tFOX_nearby\tKLF_nearby\tCTCF_nearby\tNonRespNear\tnumPeaks")
multipeaks=getMultipeaks(enhancers)
singletons=getSingletons(enhancers)
for enhancer in enhancers:
#for enhancer in multipeaks:
#for enhancer in singletons:
    extractFeatures(enhancer)


