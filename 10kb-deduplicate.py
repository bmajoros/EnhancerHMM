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

MARGIN=5000
TIMEPOINTS=("t00","t05","t1","t2","t3","t4","t5","t6","t7","t8","t10","t12")

def loadP300(filename):
    byChr={}
    with open(filename,"rt") as IN:
        for line in IN:
            fields=line.rstrip().split()
            if(len(fields)!=4): continue
            (chr,begin,end,peakID)=fields
            begin=int(begin); end=int(end)
            begin=end=int((begin+end)/2)
            array=byChr.get(chr,None)
            if(array is None): array=byChr[chr]=[]
            rec=Interval(begin,end)
            rec.chr=chr
            rec.id=peakID
            rec.link=None
            rec.linked=False
            array.append(rec)
    return byChr

def deduplicate(array,inDir,outDir):
    clusters=getClusters(array)
    reps=pickReps(clusters)
    for rep in reps:
        for time in TIMEPOINTS:
            cmd="ln -s "+inDir+"/"+rep.id+\
                ".standardized_across_all_timepoints."+time+".fastb ."
            print(cmd)

def clusterToSet(cluster):
    s=set()
    s.add(cluster)
    while(cluster.link is not None):
        s.add(cluster.link)
        cluster=cluster.link
    return s

def getMeanDist(elem,cluster):
    sum=0
    N=len(cluster)
    for other in cluster:
        sum+=elem.distance(other)
    return float(sum)/float(N)

def pickReps(clusters):
    reps=set()
    for cluster in clusters:
        asSet=clusterToSet(cluster)
        minDist=None
        rep=None
        for elem in asSet:
            dist=getMeanDist(elem,asSet)
            if(minDist is None or dist<minDist):
                minDist=dist
                rep=elem
        reps.add(rep)
    return reps

def overlap(this,next):
    return this.distance(next)<MARGIN

def getClusters(array):
    N=len(array)
    i=0
    while(i+1<N):
        this=array[i]; next=array[i+1]
        if(overlap(this,next)):
            this.link=next
            next.linked=True
        i+=1
    clusters=[]
    for i in range(N):
        elem=array[i]
        if(not elem.linked): clusters.append(elem)
    return clusters

#=========================================================================
# main()
#=========================================================================
if(len(sys.argv)!=4):
    exit(ProgramName.get()+" <p300.bed> <in-dir> <out-dir>\n")
(p300File,inDir,outDir)=sys.argv[1:]

byChr=loadP300(p300File)
chroms=byChr.keys()
for chr in chroms:
    array=byChr[chr]
    deduplicate(array,inDir,outDir)


