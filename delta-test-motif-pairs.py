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
import commands

def fisher(a,b,c,d):
    cmd="/home/bmajoros/src/scripts/fisher-exact-test.R "+str(a)+" "+\
        str(b)+" "+str(c)+" "+str(d)
    (ret,line)=commands.getstatusoutput(cmd)
    return float(line.rstrip())

#=========================================================================
# main()
#=========================================================================
if(len(sys.argv)!=2):
    exit(ProgramName.get()+" <motif-pairs.txt>\n")
infile=sys.argv[1]

peaks=[]
pairs=[]
motifs=set()
IN=open(infile,"rt")
IN.readline()
for line in IN:
    fields=line.rstrip().split()
    if(len(fields)!=2): continue
    (first,second)=fields
    if(first=="peak"):
        peaks.append(pairs)
        pairs=[]
    else: 
        pairs.append([first,second])
        motifs.add(first); motifs.add(second)
IN.close()

for motif1 in motifs:
    for motif2 in motifs:
        if(motif1>motif2): continue
        
