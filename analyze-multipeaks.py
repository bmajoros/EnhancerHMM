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
from Shuffler import Shuffler
from SummaryStats import SummaryStats
from Rex import Rex
rex=Rex()

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
            #if(peak.open and not peak.active): return True
            if(not peak.active): return True
        return False
    def hasActive(self):
        for peak in self.peaks:
            if(peak.active): return True
        return False
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
            if(not rex.find("(\d+):(\d+)-(\d)",parseElem)):
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
            peak.nearInactive=False
            peak.nearActive=False
            if(peak.active and not peak.open): activeButClosed+=1
            if(peak.active): active+=1
            peak.motifs=set()
            if(GR[i]=="1"): peak.motifs.add("GR")
            if(AP1[i]=="1"): peak.motifs.add("AP1")
            if(CEBP[i]=="1"): peak.motifs.add("CEBP")
            if(FOX[i]=="1"): peak.motifs.add("FOX")
            if(KLF[i]=="1"): peak.motifs.add("KLF")
            if(CTCF[i]=="1"): peak.motifs.add("CTCF")
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
    IN.close()
    #print("active:",active)
    #print("active but closed:",activeButClosed)
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

def analyzeInactives(multipeaks,verbose):
    withInactive=0
    withActive=0
    withBoth=0
    N=len(multipeaks)
    totalActive=0
    totalInactive=0
    for enhancer in multipeaks:
       if(enhancer.hasInactive()): withInactive+=1
       if(enhancer.hasActive()): withActive+=1
       if(enhancer.hasInactive() and enhancer.hasActive()): withBoth+=1
       #if(enhancer.hasOpen()): N+=1
       totalActive+=enhancer.numActivePeaks()
       totalInactive+=enhancer.numInactivePeaks()
    if(verbose):
        print("withInactive="+str(withInactive)+"\twithActive="+str(withActive))
        print("N="+str(N)+"\twithBoth="+str(withBoth))
        print("total active=",totalActive,"total inactive=",totalInactive)
    return (withActive,withInactive,withBoth,N)
    
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

def computePvalues(withActive,withInactive,withBoth):
    randomWithBoth=[]; randomWithInactive=[]; randomWithActive=[]
    BOTH=open("random.both","wt")
    ACTIVE=open("random.active","wt")
    INACTIVE=open("random.inactive","wt")
    for i in range(1000):
        shuffle(multipeaks)
        (withAct,withInact,withB,N)=analyzeInactives(multipeaks,False)
        randomWithBoth.append(withB)
        randomWithInactive.append(withInact)
        randomWithActive.append(withAct)
        print(withB,file=BOTH,flush=True)
        print(withAct,file=ACTIVE,flush=True)
        print(withInact,file=INACTIVE,flush=True)
    BOTH.close(); ACTIVE.close(); INACTIVE.close()
    (P,count)=getPvalue(randomWithBoth,withBoth)
    print("P-value for both:",P,str(count)+"/"+str(len(randomWithBoth)))
    (P,count)=getPvalue(randomWithActive,withActive)
    print("P-value for active:",P,str(count)+"/"+str(len(randomWithActive)))
    (P,count)=getPvalue(randomWithInactive,withInactive)
    print("P-value for inactive:",P,str(count)+"/"+str(len(randomWithInactive)))

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
        results.append(enhancer)
    return results

def getP300nearInactive(multipeaks):
    values=[]
    for enhancer in multipeaks:
        if(not enhancer.hasInactive()): continue
        for peak in enhancer.peaks:
            if(peak.active): values.append(peak.p300_3)
    return values

def getP300nearInactive(multipeaks):
    values=[]
    for enhancer in multipeaks:
        for peak in enhancer.peaks:
            if(peak.active and peak.nearInactive): values.append(peak.p300_3)
    return values

def getP300notNearInactive(multipeaks):
    values=[]
    for enhancer in multipeaks:
        for peak in enhancer.peaks:
            if(peak.active and not peak.nearInactive):
                values.append(peak.p300_3)
    return values

def getP300nearActive(multipeaks):
    values=[]
    for enhancer in multipeaks:
        for peak in enhancer.peaks:
            if(peak.active and peak.nearActive): values.append(peak.p300_3)
    return values

def getP300notNearActive(multipeaks):
    values=[]
    for enhancer in multipeaks:
        for peak in enhancer.peaks:
            if(peak.active and not peak.nearActive):
                values.append(peak.p300_3)
    return values

def writeValues(values,filename):
    OUT=open(filename,"wt")
    for value in values:
        OUT.write(str(value)+"\n")
    OUT.close()

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

# Analyze patterns
multipeaks=getMultipeaks(enhancers)
#(withActive,withInactive,withBoth,N)=analyzeInactives(multipeaks,True)
#singletons=getSingletons(enhancers)
p300Inactive=getP300nearInactive(multipeaks)
p300NoInactive=getP300notNearInactive(multipeaks)
(U,P)=mannwhitneyu(p300Inactive,p300NoInactive)
(meanInactive,SD,min,max)=SummaryStats.roundedSummaryStats(p300Inactive)
(meanNoInactive,SD,min,max)=SummaryStats.roundedSummaryStats(p300NoInactive)
print("near inactive=",meanInactive,"not near inactive=",meanNoInactive,
      "U=",U,"P=",P)
writeValues(p300Inactive,"p300.inactive")
writeValues(p300NoInactive,"p300.noinactive")

p300Active=getP300nearActive(multipeaks)
p300NoActive=getP300notNearActive(multipeaks)
(U,P)=mannwhitneyu(p300Active,p300NoActive)
(meanActive,SD,min,max)=SummaryStats.roundedSummaryStats(p300Active)
(meanNoActive,SD,min,max)=SummaryStats.roundedSummaryStats(p300NoActive)
print("near active=",meanActive,"not near active=",meanNoActive,
      "U=",U,"P=",P)
writeValues(p300Inactive,"p300.active")
writeValues(p300NoInactive,"p300.noactive")

# Shuffle to generate empirical P-values
if(WANT_SHUFFLE):
    computePvalues(withActive,withInactive,withBoth)




