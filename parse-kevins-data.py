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
import os
from Pipe import Pipe
from Rex import Rex
rex=Rex()

ALPHA=0.05
INDIR="/data/reddylab/projects/GGR/results/occupancy_prediction/hg38/jasparfix_memedb/predictions_occupancy_peaks/OpenChromDnase/1e-5/"
OUTDIR="/home/bmajoros/GGR/delta/kevin"

def parse(filename,factor):
    pipe=Pipe("cat "+filename+" | gunzip")
    pipe.readline() # header
    while(True):
        line=pipe.readline()
        if(not line): break
        fields=line.rstrip().split()
        if(len(fields)!=18): continue
        (chr,start,end,name,pwm_score,strand,t00,t05,t1,t2,t3,t4,t5,t6,
         t7,t8,t10,t12)=fields
        keep=False; score=None
        for field in fields[6:]:
            if(float(field)<=ALPHA): keep=True; score=float(field)
        if(keep):
            print(chr,start,end,factor,score,strand,sep="\t")

#=========================================================================
# main()
#=========================================================================
listing=os.listdir(INDIR)
for dir in listing:
    if(not rex.find("(\S+)_(\S+)_1e-5",dir)): continue
    factor=rex[1]; jasparID=rex[2]
    if(rex.find("acceptor",factor)): continue
    files=os.listdir(INDIR+dir)
    rankedFiles=[]
    for file in files:
        if(not rex.find("p.adjust.txt.gz",file)): continue
        level=None
        if(rex.find("Top",file)): level=1
        elif(rex.find("Middle",file)): level=2
        elif(rex.find("Bottom",file)): level=3
        else: raise Exception(file)
        rankedFiles.append({"file":file,"rank":level})
    rankedFiles.sort(key=lambda x: x["rank"])
    file=rankedFiles[0]
    parse(INDIR+dir+"/"+file["file"],factor)
