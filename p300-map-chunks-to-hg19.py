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

def readP300(filename,halfWidth):
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

#=========================================================================
# main()
#=========================================================================
if(len(sys.argv)!=4):
    exit(ProgramName.get()+" <genomewide.txt> <p300.bed> <halfWidth>\n")
(infile,p300bed,halfWidth)=sys.argv[1:]
halfWidth=int(halfWidth)

p300hash=readP300(p300bed,halfWidth)
IN=open(infile,"rt")
while(True):
    line=IN.readline()
    if(not line): break
    fields=line.split()
    if(len(fields)!=6): continue
    (name,begin,end,elemID,LLR,parse)=fields
    if(not rex.find("^(\S+).standardized",name)):
        exit("can't parse name: "+name)
    id=rex[1]
    if(p300hash.get(id,None) is None): continue #exit("can't find "+id)
    rec=p300hash[id]
    print(rec.chr,rec.interval.begin,rec.interval.end,name,LLR,sep="\t")
IN.close()



