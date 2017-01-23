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
import sys
import ProgramName
from BedReader import BedReader
from Bed3Record import Bed3Record

def readCoords(filename,halfWidth):
    hash={}
    reader=BedReader(filename)
    while(True):
        record=reader.nextRecord() # Bed3Record or Bed6Record
        if(not record): break
        center=int((record.interval.begin+record.interval.end)/2)
        rec=Bed3Record(record.chr,center-halfWidth,center+halfWidth)
        hash[record.name]=rec
    reader.close()
    return hash

def inc(key,hash):
    hash[key]=hash.get(key,0)+1

def dumpHash(hash,label):
    print(label+":")
    keys=hash.keys()
    for key in keys:
        print("\t"+key+" = "+str(hash[key]))

def printElem(halfHour,oneHour,coord,elemID):
    type=str(halfHour)+str(oneHour)
    print(coord.toString()+"\t"+elemID+"_"+type)

#=========================================================================
# main()
#=========================================================================
if(len(sys.argv)!=3):
    exit(ProgramName.get()+" <trajectories-num-peaks.txt> <p300.bed>\n")
(trajectoriesFile,coordsFile)=sys.argv[1:]

coords=readCoords(coordsFile,1000)
IN=open(trajectoriesFile,"rt")
changeAt1={}; changeAt2={}
for line in IN:
    fields=line.rstrip().split()
    if(len(fields)!=12): continue
    (id,halfHour,oneHour,twoHour)=fields[:4]
    inc(halfHour+"->"+oneHour,changeAt1)
    inc(halfHour+"->"+twoHour,changeAt2)
    halfHour=int(halfHour); oneHour=int(oneHour); twoHour=int(twoHour)
    coord=coords[id]
    if(halfHour==2 and oneHour==0): printElem(2,0,coord,id)
    elif(halfHour==0 and oneHour==1): printElem(0,1,coord,id)
IN.close()
#dumpHash(changeAt1,"one hour")
#dumpHash(changeAt2,"two hour")



