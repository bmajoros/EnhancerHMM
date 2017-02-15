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

def loadLoopFile(filename):
    # ENSG00000164880.14      iter0_peak16234
    with open(filename) as IN:
        for line in IN:
            fields=line.rstrip().split()
            if(len(fields)!=2): continue
            (gene,peak)=fields
            

def loadP300file(filename):
    # iter0_peak14413.t3.fastb        2234.2906200000034      1.0     1:0-1|2:1-2|3:2-1243|4:1243-1256|3:1256-1998|4:1998-1999|5:1999-2000    1.07942122577,1.12211048111,1.00477844983,1.16141454506|0.952633976593,1.02660980228,0.892990209161,0.961271510352      10      10      11      10      10      00
    pass

def loadExpressionFile(filename):
    # gene        logFC   logCPM  LR      PValue  FDR
    pass

#=========================================================================
# main()
#=========================================================================
if(len(sys.argv)!=4):
    exit(ProgramName.get()+" <genes-enhancers.txt> <genomewide-features.txt> <edgeR.txt>\n")
(loopFile,p300file,expressionFile)=sys.argv[1:]






