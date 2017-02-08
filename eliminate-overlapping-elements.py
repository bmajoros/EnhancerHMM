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
from BedReader import BedReader
import os
from Rex import Rex
rex=Rex()

DELTA="/home/bmajoros/GGR/delta"
P300_BED=DELTA+"/p300-hg19.bed"
PARTITIONS=DELTA+"/fastb/crossval/partitions" # /pos/1
TRAINING_SETS=DELTA+"/fastb/crossval/training-sets" # /pos/leaveout1

def getBlacklist(p300):
    blacklist=set()
    chroms=p300.keys()
    for chr in chroms:
        array=p300[chr]
        for rec in array:
            interval=rec.interval
            center=interval.intCenter()
            interval.begin=center-1000
            interval.end=center+1000
        array.sort(key=lambda x: x.interval.begin)
        N=len(array)
        for i in range(N-1):
            rec=array[i]
            nextRec=array[i+1]
            if(rec.interval.overlaps(nextRec.interval)):
                blacklist.add(rec.name)
    return blacklist

def filter(dir,subs,blacklist):
    for sub in subs:
        files=os.listdir(dir+"/"+sub)
        for file in files:
            if(not rex.find("([^\/]+)\.fastb",file)): continue
            id=rex[1]
            if(id in blacklist):
                cmd="mv "+dir+"/"+sub+"/"+file+" "+dir+"/removed"
                os.system(cmd)

#=========================================================================
# main()
#=========================================================================

p300=BedReader.hashBySubstrate(P300_BED) # chr -> list of records
blacklist=getBlacklist(p300)
filter(PARTITIONS+"/pos",["1","2","3","4","5"],blacklist)
filter(PARTITIONS+"/neg",["1","2","3","4","5"],blacklist)
filter(TRAINING_SETS+"/pos",["leaveout1","leaveout2","leaveout3","leaveout4",
                             "leaveout5"],blacklist)
filter(TRAINING_SETS+"/neg",["leaveout1","leaveout2","leaveout3","leaveout4",
                             "leaveout5"],blacklist)





