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
from Interval import Interval
rex=Rex()

def loadPeaks(filename):
    array=[]
    with open(filename,"rt") as IN:
        for line in IN:
            fields=line.split()
            if(len(fields)<10): continue
            if(not rex.find("_peak_",fields[9])): continue
            (chr,start,end,length,abs_summit,pileup,pvalue,fold_enrichment,
             qvalue,name)=fields
            elem=Interval(int(start),int(end))
            elem.chr=chr
            elem.pvalue=pvalue
            elem.qvalue=qvalue
            elem.fold=fold_enrichment
            array.append(elem)
    array.sort(key=lambda i: i.begin)
    return array

#=========================================================================
# main()
#=========================================================================
if(len(sys.argv)!=3):
    exit(sys.argv[0]+" <dex.xls> <eth.xls>")
(dexFile,ethFile)=sys.argv[1:]

# Load peaks
dexPeaks=loadPeaks(dexFile)
ethPeaks=loadPeaks(ethFile)

for dexPeak in dexPeaks:
    overlaps=False
    for ethPeak in ethPeaks:
        if(ethPeak.overlaps(dexPeak)): 
            overlaps=True
            break
    if(not overlaps):
        print(dexPeak.chr,dexPeak.begin,dexPeak.end,"\t")




