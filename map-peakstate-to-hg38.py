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
from Rex import Rex
rex=Rex()

def readP300(filename):
    hash={}
    reader=BedReader(filename)
    while(True):
        record=reader.nextRecord() # Bed3Record or Bed6Record
        if(not record): break
        hash[record.name]=record
    reader.close()
    return hash

#=========================================================================
# main()
#=========================================================================
if(len(sys.argv)!=5):
    exit(ProgramName.get()+" <genomewide.txt> <p300.bed> <posterior-threshold> <min-len>\n")
(infile,p300bed,threshold,MIN_LEN)=sys.argv[1:]
THRESHOLD=float(threshold)
MIN_LEN=int(MIN_LEN)

p300hash=readP300(p300bed)
IN=open(infile,"rt")
while(True):
    line=IN.readline()
    if(not line): break
    fields=line.split()
    if(len(fields)!=8): continue
    (name,begin,end,elemID,LLR,posterior,parse,topBottom)=fields
    if(float(posterior)<THRESHOLD): continue
    id=None
    if(rex.find("^(\S+).standardized",name)): id=rex[1]
    elif(rex.find("^(\S+)\.t\d+\.fastb",name)): id=rex[1]
    else: exit("can't parse name: "+name)
    if(p300hash.get(id,None) is None): continue #exit("can't find "+id)
    rec=p300hash[id]
    origin=rec.interval.begin-1000
    subfields=parse.split("|")
    for field in subfields:
        if(not rex.find("(\d+):(\d+)-(\d+)",field)):
            raise Exception("can't parse "+field)
        state=int(rex[1])
        if(state!=3): continue
        begin=int(rex[2]); end=int(rex[3]); L=end-begin
        if(L<MIN_LEN): continue
        print(rec.chr,origin+begin,origin+end,sep="\t")
IN.close()



