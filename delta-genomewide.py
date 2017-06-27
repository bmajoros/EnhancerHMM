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
from Fastb import Fastb
from FastbTrack import FastbTrack
from Rex import Rex
rex=Rex()

# GLOBALS
BASE="/home/bmajoros/GGR/delta"
#NODNA=BASE+"/no-dna"
NODNA=BASE+"/10kb-raw"
#OUTDIR=BASE+"/delta-fastb"
OUTDIR=BASE+"/10kb"
TIMEPOINTS=("t05","t1","t2","t3","t4","t5","t6","t7","t8","t10","t12")

def loadList(infile):
    peaks=[]
    with open(infile,"rt") as IN:
        for line in IN:
            line=line.strip()
            if(line!=""):
                peaks.append(line)
    return peaks

def process(peaks,indir,outdir):
    for peak in peaks:
        for timepoint in TIMEPOINTS:
            in1=indir+"/"+peak+".standardized_across_all_timepoints.t00.fastb"
            in2=indir+"/"+peak+".standardized_across_all_timepoints."+\
                timepoint+".fastb"
            if(not os.path.exists(in1) or not os.path.exists(in2)): continue
            outfile=outdir+"/"+peak+"."+timepoint+".fastb"
            combine(in1,in2,timepoint,outfile)

def combine(in1,in2,timepoint,outfile):
    fastb1=Fastb(in1);
    fastb2=Fastb(in2);
    n=fastb1.numTracks()
    toDrop=set()
    for i in range(n):
        track=fastb1.getIthTrack(i)
        name=track.getID()
        if(track.isDiscrete()): toDrop.add(name)
        else: fastb1.renameTrack(name,name+".t00")
    for d in toDrop: fastb1.dropTrack(d)
    n=fastb2.numTracks()
    for i in range(n):
        track=fastb2.getIthTrack(i)
        name=track.getID()
        if(track.isDiscrete()): continue
        track.rename(name+".t3") # name used by the HMM
        fastb1.addTrack(track);
    fastb1.save(outfile)

#=========================================================================
# main()
#=========================================================================
if(len(sys.argv)!=2):
    exit(ProgramName.get()+" <partition.txt>\n")
infile=sys.argv[1]
peaks=loadList(infile)
process(peaks,NODNA,OUTDIR)


