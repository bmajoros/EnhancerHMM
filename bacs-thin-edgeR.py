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
from Rex import Rex
rex=Rex()

def updateHash(line,fields,hash):
    if(not rex.find("(\S+)_(\d+)$",fields[0])):
        raise Exception("Can't parse first column: "+fields[0])
    chr=rex[1]
    pos=int(rex[2])
    if(not hash.get(chr,None)): hash[chr]=[]
    array=hash[chr]
    array.append(line)
    #fields=list(map(lambda x:int(x),fields[1:]))
    #array.append(fields)
    #if(pos!=len(array)): exit(str(pos)+"!="+str(len(array)))

#=========================================================================
#                                main()
#=========================================================================
if(len(sys.argv)!=4):
    exit(sys.argv[0]+" <in.txt> <distance> <out.txt>")
(infile,distance,outfile)=sys.argv[1:]
distance=int(distance)

# Read table
hash={}
IN=open(infile,"rt")
header=IN.readline()
while(True):
    line=IN.readline()
    if(not line): break
    fields=line.split()
    if(len(fields)!=8): continue
    updateHash(line,fields,hash)
IN.close()

OUT=open(outfile,"wt")
OUT.write(header)
keys=hash.keys()
for chr in keys:
    array=hash[chr]
    L=len(array)
    for i in range(0,L,distance):
        line=array[i]
        OUT.write(line)
OUT.close()
    

