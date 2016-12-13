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

# Process command line
if(len(sys.argv)!=4):
    exit(sys.argv[0]+" <in.path> <substrate> <out.bed>\n")
(infile,chr,outfile)=sys.argv[1:]

# Read path file
IN=open(infile,"rt")
header=IN.readline()
array=[]
while(True):
    line=IN.readline()
    if(not line): break
    state=int(line)
    array.append(state)
IN.close()

# Find predicted elements
OUT=open(outfile,"wt")
L=len(array)
begin=-1
for i in range(L):
    if(array[i]>1):
        if(i>0 and array[i-1]==1 or i==0):
            begin=i
    else:
        if(i>0 and array[i-1]>1):
            OUT.write(chr+"\t"+str(begin)+"\t"+str(i)+"\n")
OUT.close()


