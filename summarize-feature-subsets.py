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
from SummaryStats import SummaryStats

KEYS=("full_model",
      "all_t00",
      "all_t3",
      "H3K2AC_both_times",
      "P300_both_times",
      "H3K4ME1_both_times",
      "H3K4ME2_both_times",
      "DNASE_both_times",
      "DNase.t00",
      "H3K27ac.t00",
      "H3K4me1.t00",
      "H3K4me2.t00",
      "EP300.t00",
      "DNase.t3",
      "H3K27ac.t3",
      "H3K4me1.t3",
      "H3K4me2.t3",
      "EP300.t3"
      )

#=========================================================================
# main()
#=========================================================================
if(len(sys.argv)!=2):
    exit(ProgramName.get()+" <test-feature-subsets.txt>\n")
infile=sys.argv[1]

hash={}
with open(infile,"rt") as IN:
    for line in IN:
        fields=line.rstrip().split()
        if(len(fields)!=3): continue
        (EVAL,key,score)=fields
        if(EVAL!="EVAL"): continue
        if(hash.get(key,None) is None): hash[key]=[]
        hash[key].append(float(score))
print("features\tmean\tSD\tmin\tmax")
for key in KEYS:
    [mean,SD,min,max]=SummaryStats.roundedSummaryStats(hash[key])
    print(key,mean,SD,min,max,sep="\t")




