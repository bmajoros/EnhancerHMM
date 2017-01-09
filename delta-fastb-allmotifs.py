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
from Rex import Rex
rex=Rex()

# GLOBALS
BASE="/home/bmajoros/GGR/delta"
TEST_POS_FROM=BASE+"/test-pos-nomotif"
TEST_NEG_FROM=BASE+"/test-neg-nomotif"
TEST_POS_TO=BASE+"/test-pos-allmotifs"
TEST_NEG_TO=BASE+"/test-neg-allmotifs"
TRAIN_POS_FROM=BASE+"/train-pos"
TRAIN_NEG_FROM=BASE+"/train-neg"
TRAIN_POS_TO=BASE+"/train-pos-allmotifs"
TRAIN_NEG_TO=BASE+"/train-neg-allmotifs"
MOTIF_TRACKS=BASE+"/motifs-discrete"
MUMMIE=os.environ["MUMMIE"]
TRACKS="DNase.t00 H3K27ac.t00 H3K4me1.t00 H3K4me2.t00 EP300.t00 DNase.t3 H3K27ac.t3 H3K4me1.t3 H3K4me2.t3 EP300.t3"

def processTest(fromDir,toDir):
    files=os.listdir(fromDir)
    for file in files:
        continuousFile=fromDir+"/"+file
        if(not rex.find("(\S+).fastb",file)):
            exit("can't parse filename "+file)
        filestem=rex[1]
        toFile=toDir+"/"+file
        motifFile=MOTIF_TRACKS+"/"+filestem+\
            ".standardized_across_all_timepoints.t00.fastb"
        cmd="cat "+continuousFile+" "+motifFile+">"+toFile
        os.system(cmd)

def processTrain(fromDir,toDir):
    files=os.listdir(fromDir)
    for file in files:
        continuousFile=fromDir+"/"+file
        if(not rex.find("(\S+).fastb",file)):
            exit("can't parse filename "+file)
        filestem=rex[1]
        toFile=toDir+"/"+file
        motifFile=MOTIF_TRACKS+"/"+filestem+\
            ".standardized_across_all_timepoints.t00.fastb"
        os.system(MUMMIE+"/fastb-extract-tracks.pl "+continuousFile+
                  " tmp.1 "+TRACKS)
        cmd="cat tmp.1 "+motifFile+">"+toFile
        os.system(cmd)

#=========================================================================
# main()
#=========================================================================
#processTest(TEST_POS_FROM,TEST_POS_TO)
#processTest(TEST_NEG_FROM,TEST_NEG_TO)
processTrain(TRAIN_POS_FROM,TRAIN_POS_TO)
processTrain(TRAIN_NEG_FROM,TRAIN_NEG_TO)

