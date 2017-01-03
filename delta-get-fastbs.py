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
import TempFilename

# GLOBALS
BASE="/home/bmajoros/GGR/delta"
POS_LIST=BASE+"/positives.txt"
NEG_LIST=BASE+"/negatives.txt"
POS_DIR=BASE+"/positives"
NEG_DIR=BASE+"/negatives"
NODNA=BASE+"/no-dna"
MOTIF_TRACKS=BASE+"/motif-tracks"
MUMMIE=os.environ["MUMMIE"]
tempFile1=TempFilename.generate(".fastb")
tempFile2=TempFilename.generate(".fastb")
tempFile3=TempFilename.generate(".fastb")

def loadList(filename):
    array=[]
    IN=open(filename,"rt")
    for line in IN:
        array.append(line.rstrip())
    IN.close()
    return array

def process(listfile,outdir):
    files=loadList(listfile)
    for file in files:
        copyTo(file,outdir)

def getMotifTrack(peak,outfile):
    program=MUMMIE+"/fastb-extract-tracks.pl"
    motifFastb=MOTIF_TRACKS+"/"+peak+\
        ".standardized_across_all_timepoints.t00.fastb"
    cmd=program+" "+motifFastb+" "+outfile+" GR/AR/MR"
    os.system(cmd)

def renameTrack(filename,oldName,newName):
    cmd=MUMMIE+"/fastb-rename-track.pl "+filename+" "+oldName+" "+newName
    os.system(cmd)

def getBaseTracks(peak,timepoint,outfile):
    baseFastb=NODNA+"/"+peak+".standardized_across_all_timepoints."+timepoint\
        +".fastb"
    os.system("cp "+baseFastb+" "+outfile)
    renameTrack(outfile,"DNase","DNase."+timepoint)
    renameTrack(outfile,"EP300","EP300."+timepoint)
    renameTrack(outfile,"H3K27ac","H3K27ac."+timepoint)
    renameTrack(outfile,"H3K4me1","H3K4me1."+timepoint)
    renameTrack(outfile,"H3K4me2","H3K4me2."+timepoint)

def copyTo(peak,outdir):
    getMotifTrack(peak,tempFile1)
    getBaseTracks(peak,"t00",tempFile2)
    getBaseTracks(peak,"t3",tempFile3)
    cmd="cat "+tempFile2+" "+tempFile3+" "+tempFile1+" > "+outdir+"/"+peak\
        +".fastb"
    os.system(cmd)

#=========================================================================
# main()
#=========================================================================
process(POS_LIST,POS_DIR)
process(NEG_LIST,NEG_DIR)
os.remove(tempFile1)
os.remove(tempFile2)
os.remove(tempFile3)


