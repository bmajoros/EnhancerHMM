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
from scipy.stats import fisher_exact
from Rex import Rex
rex=Rex()

ALL_MOTIFS=set()
MOTIFS=("CTCF","AP1","GR","FOX","KLF","CEBP")
for motif in MOTIFS: ALL_MOTIFS.add(motif)

#=========================================================================
# main()
#=========================================================================
if(len(sys.argv)!=2):
    exit(ProgramName.get()+" <infile>\n")
(infile,)=sys.argv[1:]

peakHits={}; peakMisses={}; humpHits={}; humpMisses={}
with open(infile,"rt") as IN:
    for line in IN:
        fields=line.rstrip().split()
        if(len(fields)!=9): raise Exception(line)
        peakHump=fields[1]
        hits=set()
        for pair in fields[3:]:
            if(not rex.find("(\S+)=(\d+)",pair)): raise Exception(pair)
            motif=rex[1]; count=int(rex[2])
            if(count>0): hits.add(motif)
        for motif in ALL_MOTIFS:
            if(motif in hits):
                if(peakHump=="peak"): peakHits[motif]=peakHits.get(motif,0)+1
                else: humpHits[motif]=humpHits.get(motif,0)+1
            else:
                if(peakHump=="peak"): 
                    peakMisses[motif]=peakMisses.get(motif,0)+1
                else: humpMisses[motif]=humpMisses.get(motif,0)+1
for motif in ALL_MOTIFS:
    peakCount=peakHits[motif]; peakAbsent=peakMisses[motif]
    humpCount=humpHits[motif]; humpAbsent=humpMisses[motif]
    (oddsRatio,P)=fisher_exact([[peakCount,peakAbsent],[humpCount,humpAbsent]])
    enrichment=(float(peakCount)/float(peakCount+peakAbsent))/\
        (float(humpCount)/float(humpCount+humpAbsent))
    print(enrichment,P,motif,peakCount,peakAbsent,humpCount,humpAbsent,sep="\t")

# iter0_peak10009_1       hump    100     AP1=0   CTCF=0  FOX=0   CEBP=0  KLF=0   GR=0


