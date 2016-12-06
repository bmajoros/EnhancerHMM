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

# Globals
BASE="/home/bmajoros/GGR/STAR/output"
DEPTH_FILES=("dna-depth.bed",
             "3hr-EtOH-Rep1-depth.bed",
             "3hr-EtOH-Rep2-depth.bed",
             "3hr-EtOH-Rep3-depth.bed",
             "3hr-Dex-Rep1-depth.bed",
             "3hr-Dex-Rep2-depth.bed",
             "3hr-Dex-Rep3-depth.bed")
OUT_FILE="edgeR-input.txt"

# Functions
def nextLine(handles,OUT):
    label=None
    outline=""
    for handle in handles:
        line=handle.readline()
        if(not line): return None
        fields=line.split()
        if(len(fields)!=3): raise Exception("can't parseline: "+line)
        (chr,pos,reads)=fields
        if(not label): label=chr+"_"+pos
        outline+="\t"+reads
    OUT.write(label+outline+"\n")
    return True

#=========================================================================
#                                 main()
#=========================================================================
OUT=open(OUT_FILE,"wt")
OUT.write("dna\tEtOH1\tEtOH2\tEtOH3\tDex1\tDex2\tDex3\n")
handles=[]
for file in DEPTH_FILES:
    handles.append(open(file,"rt"))
while(True):
    if(not nextLine(handles,OUT)): break
for handle in handles: handle.close()
OUT.close()

