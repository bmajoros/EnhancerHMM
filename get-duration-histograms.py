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
import ProgramName
from Rex import Rex
rex=Rex()

#=========================================================================
# main()
#=========================================================================
if(len(sys.argv)!=3):
    exit(ProgramName.get()+" <genomewide.txt> <LLR-threshold>\n")
(predictionsFile,threshold,)=sys.argv[1:]
threshold=float(threshold)

files=[None]*5
for i in range(5):
    files[i]=open("state"+str(i+1)+".durations","wt")
with open(predictionsFile) as IN:
    for line in IN:
        fields=line.rstrip().split()
        if(len(fields)!=8): continue
        (fastb,b,e,taskID,LLR,P,parse,topBottom)=fields
        if(float(LLR)<threshold): continue
        #fields=parse.split("|")
        #for field in fields:
        #    if(not rex.find("(\d+):(\d+)-(\d+)",field)): exit(field)
        #    state=int(rex[1]); begin=int(rex[2]); end=int(rex[3])
        #    L=end-begin
        #    print(L,file=files[state-1])
        if(not rex.find("(\d+)-(\d+):(\d+)-(\d+):(\d+)-(\d+)",parse)):
            exit(parse)
        print(int(rex[1]),file=files[0])
        print(int(rex[2])-int(rex[1]),file=files[1])
        print(int(rex[4])-int(rex[3]),file=files[2])
        print(int(rex[6])-int(rex[5]),file=files[3])
        print(2000-int(rex[6]),file=files[4])
for file in files: file.close()



