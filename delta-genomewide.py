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
NODNA=BASE+"/no-dna"
OUTDIR=BASE+"/delta-fastb"
TIMEPOINTS=("t05","t1","t2","t3","t4","t5","t6","t7","t8","t10","t12")

def loadList(infile):
    peaks=[]
    with open(infile,"rt") as IN:
        for line in IN:
            line=line.strip()
            if(line!=""):
                peaks.append(line)
    return peaks

#def getPeakList(dir):
#    peaks=[]
#    files=os.listdir(dir)
#    for file in files:
#        if(rex.find("(\S+)\.standard.*\.t00",file)):
#            peaks.append(rex[1])
#    return peaks

def process(peaks,indir,outdir):
    #peaks=getPeakList(indir)
    for peak in peaks:
        for timepoint in TIMEPOINTS:
            in1=indir+"/"+peak+".standardized_across_all_timepoints.t00.fastb"
            in2=indir+"/"+peak+".standardized_across_all_timepoints."+\
                timepoint+".fastb"
            outfile=outdir+"/"+peak+"."+timepoint+".fastb"
            combine(in1,in2,timepoint,outfile)

def combine(in1,in2,timepoint,outfile):
    fastb1=Fastb(in1);
    fastb2=Fastb(in2);
    n=fastb1.numTracks()
    for i in range(n):
        track=fastb1.getIthTrack(i)
        name=track.getID()
        fastb1.renameTrack(name,name+".t00")
    n=fastb2.numTracks()
    for i in range(n):
        track=fastb2.getIthTrack(i)
        name=track.getID()
        track.rename(name+"."+timepoint)
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


