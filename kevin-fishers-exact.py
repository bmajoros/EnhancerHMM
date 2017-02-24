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
from Fisher2x2 import Fisher2x2

def loadData(filename):
    lines=[]
    IN=open(filename)
    for line in IN:
        fields=line.rstrip().split()
        rec=[fields[0],set()]
        for TF in fields[1:]: rec[1].add(TF)
        lines.append(rec)
    IN.close()
    return lines

def getFactors(data):
    factors=set()
    for line in data:
        for TF in line[1]: factors.add(TF)
    return factors

def testFactors(factors,data,alpha):
    for TF in factors:
        testFactor(TF,data,alpha)

def testFactor(TF,data,alpha):
    table=[[0,0],[0,0]]
    for line in data:
        (peakType,factors)=line
        peakType=0 if peakType=="singleton" else 1
        present=1 if TF in factors else 0
        table[peakType][present]+=1
    fisher=Fisher2x2(table[0][0],table[0][1],table[1][0],table[1][1])
    (exp00,exp01,exp10,exp11)=fisher.getExpectedCounts()
    if(table[1][1]<=exp11): return
    P=fisher.getPvalue()
    if(P>=alpha): return
    print(P,"\t"+TF+"\t"+"observed=[",
          table[0][0],table[0][1],table[1][0],table[1][1],
          "] expected= [",exp00,exp01,exp10,exp11,"]")

#=========================================================================
# main()
#=========================================================================
if(len(sys.argv)!=3):
    exit(ProgramName.get()+" <motifs-in-peaks.txt> <alpha>\n")
(infile,alpha)=sys.argv[1:]
alpha=float(alpha)

data=loadData(infile)
factors=getFactors(data)
testFactors(factors,data,alpha)







