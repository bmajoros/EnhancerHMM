#!/usr/bin/env python
#=========================================================================
# This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
# License (GPL) version 3, as described at www.opensource.org.
# Copyright (C)2016 William H. Majoros (martiandna@gmail.com).
#=========================================================================
from __future__ import (absolute_import, division, print_function, 
   unicode_literals, generators, nested_scopes, with_statement)
from builtins import (bytes, dict, int, list, object, range, str, ascii,
   chr, hex, input, next, oct, open, pow, round, super, filter, map, zip)
# The above imports should allow this program to run in both Python 2 and
# Python 3.  You might need to update your version of module "future".
import os
import sys
import math
from BedReader import BedReader
from Pipe import Pipe
import SumLogProbs
from Interval import Interval
import ProgramName
import TempFilename
from Fastb import Fastb
from Rex import Rex
rex=Rex()
MUMMIE=os.environ["MUMMIE"]

MOTIFS_DIR="/home/bmajoros/GGR/delta/motifs-continuous"
RAW="/data/reddylab/projects/GGR/subprojects/hmm/data/EP300_not_standardized"
PEAK_STATE=3
fgPrior=math.log(0.5)
bgPrior=math.log(0.5)

def getLL(fastb,hmm):
    cmd=MUMMIE+"/get-likelihood "+hmm+" "+fastb
    pipe=Pipe(cmd)
    line=pipe.readline()
    return float(line)

def getPosterior(fastb):
    fgLL=getLL(fastb,fgHMM)
    bgLL=getLL(fastb,bgHMM)
    jointBG=bgLL+bgPrior
    jointFG=fgLL+fgPrior
    denom=SumLogProbs.sumLogProbs2(jointBG,jointFG)
    posterior=math.exp(jointFG-denom)
    return posterior

def getLLR(fastb):
    fgLL=getLL(fastb,fgHMM)
    bgLL=getLL(fastb,bgHMM)
    return fgLL-bgLL

def getPath(hmm,fastb):
    array=[]
    cmd=MUMMIE+"/parse "+hmm+" "+fastb
    pipe=Pipe(cmd)
    pipe.readline() # header
    while(True):
        line=pipe.readline()
        if(not line): break
        if(not rex.find("\S+",line)): continue
        state=int(line)
        array.append(state)    
    return array

def getSections(path):
    sections=[]
    L=len(path)
    elem=[]
    for i in range(L):
        state=path[i]
        eLen=len(elem)
        if(eLen==0 or elem[eLen-1]==state): elem.append(state)
        else:
            sections.append(elem)
            elem=[state]
    if(len(elem)>0): sections.append(elem)
    return sections

def process(dir,posHMM,negHMM,label):
  files=os.listdir()
  n=len(files)
  for i in range(n):
    file=files[i]
    file=file.rstrip()
    if(not rex.find(".fastb",file)): continue
    numer=getLL(file,posHMM)
    denom=getLL(file,negHMM)
    ratio=numer-denom
    ROC.write(str(ratio)+"\t"+str(label)+"\n")

def getParse(path):
    sections=getSections(path)
    parse=""
    pos=0
    for section in sections:
        L=len(section)
        type=section[0]
        if(len(parse)>0): parse+="|"
        begin=pos
        end=pos+L
        parse+=str(type)+":"+str(begin)+"-"+str(end)
        pos+=L
    return parse

def getParse_old(path):
    sections=getSections(path)
    if(len(sections)!=5):
        print(path)
        print(sys.argv)
        exit("cannot find five states in parse")
    parse=""
    pos=0
    for section in sections:
        L=len(section)
        type=section[0]
        if(type>1 and type<5 or type>6 and type<10 or type>11 and type<15):
            if(len(parse)>0): parse+=":"
            begin=pos
            end=pos+L
            parse+=str(begin)+"-"+str(end)
        pos+=L
    return parse

def getPeaks(parse):
    fields=parse.split("|")
    peaks=[]
    for field in fields:
        if(not rex.find("(\d+):(\d+)-(\d+)",field)):
            raise Exception("can't parse "+field)
        state=int(rex[1])
        if(state!=PEAK_STATE): continue
        begin=int(rex[2]); end=int(rex[3])
        peaks.append(Interval(begin,end))
    return peaks

def getFeature(fastb,interval,label):
    track=fastb.getTrackByName(label)
    return str(track.getMax(interval)[0])

def getFeatures_raw(peaks,dir,file):
    if(not rex.find("(\S+)\.(t\d+)\.fastb",file)):
        raise Excpetion(file)
    id=rex[1]; time=rex[2]
    if(not os.path.exists(dir+"/"+id+".t00.fastb") or
       not os.path.exists(dir+"/"+file)): return None
    fastb0=Fastb(dir+"/"+id+".t00.fastb")
    fastb1=Fastb(dir+"/"+file)
    features=""
    for peak in peaks:
        if(len(features)>0): features+="|"
        features+=getFeature(fastb0,peak,"DNase")+","
        features+=getFeature(fastb1,peak,"DNase")+","
        features+=getFeature(fastb0,peak,"EP300")+","
        features+=getFeature(fastb1,peak,"EP300")
    return features

def getFeatures_standardized(peaks,dir,file):
    if(not os.path.exists(dir+"/"+file)): return None
    fastb=Fastb(dir+"/"+file)
    features=""
    for peak in peaks:
        if(len(features)>0): features+="|"
        features+=getFeature(fastb,peak,"DNase.t00")+","
        features+=getFeature(fastb,peak,"DNase.t3")+","
        features+=getFeature(fastb,peak,"EP300.t00")+","
        features+=getFeature(fastb,peak,"EP300.t3")
    return features

def getMotif(file,peaks,label):
    if(not rex.find("(\S+)\.t",file)): raise Exception(file)
    id=rex[1]
    motifFile=MOTIFS_DIR+"/"+id+".standardized_across_all_timepoints.t00.fastb"
    fastb=Fastb(motifFile)
    track=fastb.getTrackByName(label)
    nonzeros=track.getNonzeroRegions()
    result=[0]*len(peaks)
    for i in range(len(peaks)):
        peak=peaks[i]
        for motif in nonzeros:
            if(motif.overlaps(peak)):
                result[i]=1
                break
    s=""
    for r in result: s+=str(r)
    return s

#=========================================================================
# main()
#=========================================================================
if(len(sys.argv)!=6):
    exit(ProgramName.get()+
         " <list.txt> <fastb-dir> <pos.hmm> <neg.hmm> <taskID>\n")
(fileList,inDir,fgHMM,bgHMM,taskID)=sys.argv[1:]
files=[]
with open(fileList,"rt") as IN:
    for line in IN:
        files.append(line.rstrip())
#nextID=1
for file in files:
    fullPath=inDir+"/"+file
    path=getPath(fgHMM,fullPath)
    parse=getParse(path)
    llr=getLLR(fullPath)
    P=getPosterior(fullPath)
    peaks=getPeaks(parse)
    #features=getFeatures_raw(peaks,RAW,file)
    features=getFeatures_standardized(peaks,inDir,file)
    if(features is None): continue
    gr=getMotif(file,peaks,"GR/AR/MR")
    ap1=getMotif(file,peaks,"AP1")
    cebp=getMotif(file,peaks,"CEBP")
    fox=getMotif(file,peaks,"FOX")
    klf=getMotif(file,peaks,"KLF")
    ctcf=getMotif(file,peaks,"CTCF")
    #id="task"+taskID+"_elem"+str(nextID)
    print(file+"\t"+str(llr)+"\t"+str(P)+"\t"+parse+"\t"+features+"\t"+
          gr+"\t"+ap1+"\t"+cebp+"\t"+fox+"\t"+klf+"\t"+ctcf,flush=True)
    #nextID+=1

#########################################################################
