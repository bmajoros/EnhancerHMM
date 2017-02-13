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
import os
import ProgramName
from SummaryStats import SummaryStats
from Rex import Rex
rex=Rex()

PREFIX="loopy-"
TIMEPOINTS=("t05","t1","t2","t3","t4","t5","t6","t7","t8","t10","t12")
TIMES=(0.5,1,2,3,4,5,6,7,8,10,12)

def getFactors(files):
    factors=set()
    for file in files:
        if(rex.find("-(\S+)-t\d+\.txt",file)):
            factor=rex[1]
            factors.add(factor)
    return factors

def getMeans(factor,time,indir):
    file=indir+"/"+PREFIX+factor+"-"+time+".txt"
    IN=open(file,"rt")
    singles=[]
    multis=[]
    for line in IN:
        fields=line.rstrip().split()
        if(len(fields)!=3): continue
        (score,peaks,peakLen)=fields
        array=singles if int(peaks)==1 else multis
        array.append(float(score))
    IN.close()
    #print("singles",len(singles),file)
    #print("multis",len(multis))
    (meanSingle,sd,min,max)=SummaryStats.summaryStats(singles)
    (meanMulti,sd,min,max)=SummaryStats.summaryStats(multis)
    return (meanSingle,meanMulti)

#=========================================================================
# main()
#=========================================================================
if(len(sys.argv)!=2):
    exit(ProgramName.get()+" <chip-dir>\n")
(indir,)=sys.argv[1:]

files=os.listdir(indir)
factors=getFactors(files)
for factor in factors:
    outfile1=factor+".singlepeak"
    outfile2=factor+".multipeak"
    OUT1=open(outfile1,"wt")
    OUT2=open(outfile2,"wt")
    for i in range(len(TIMEPOINTS)):
        timepoint=TIMEPOINTS[i]
        time=TIMES[i]
        (meanSingle,meanMulti)=getMeans(factor,timepoint,indir)
        print(time,meanSingle,sep="\t",file=OUT1)
        print(time,meanMulti,sep="\t",file=OUT2)
    OUT1.close(); OUT2.close()
    

