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
import commands
import ProgramName

def fisher(a,b,c,d):
    cmd="/home/bmajoros/src/scripts/fisher-exact-test.R "+str(a)+" "+\
        str(b)+" "+str(c)+" "+str(d)
    #print(cmd)
    (ret,line)=commands.getstatusoutput(cmd)
    return float(line.rstrip())

#=========================================================================
# main()
#=========================================================================

if(len(sys.argv)!=3):
    exit(ProgramName.get()+" <peaks-path.txt> <peaks-motifs.txt>\n")
(pathFile,motifFile)=sys.argv[1:]

paths={}
with open(pathFile,"rt") as IN:
    for line in IN:
        fields=line.rstrip().split()
        if(len(fields)!=2): continue
        (peak,path)=fields
        paths[peak]=path

motifList=set()
motifs={}
with open(motifFile,"rt") as IN:
    for line in IN:
        fields=line.rstrip().split()
        if(len(fields)!=2): continue
        (peak,motif)=fields
        if(motifs.get(peak,None) is None): motifs[peak]=set()
        motifs[peak].add(motif)
        motifList.add(motif)

peaks=paths.keys()
for motif in motifList:
    motifOnePeak=0; motifTwoPeaks=0; nomotifOnePeak=0; nomotifTwoPeaks=0
    for peak in peaks:
        path=paths[peak]
        present=motif in motifs.get(peak,set())
        if(present):
            if(path=="onepeak"): motifOnePeak+=1
            else: motifTwoPeaks+=1
        else:
            if(path=="onepeak"): nomotifOnePeak+=1
            else: nomotifTwoPeaks+=1
    P=fisher(motifOnePeak,motifTwoPeaks,nomotifOnePeak,nomotifTwoPeaks)
    total=motifOnePeak+motifTwoPeaks+nomotifOnePeak+nomotifTwoPeaks
    totalMotif=motifOnePeak+motifTwoPeaks
    totalOnePeak=motifOnePeak+nomotifOnePeak
    P_motif=float(totalMotif)/float(total)
    P_nomotif=1.0-P_motif
    P_onePeak=float(totalOnePeak)/float(total)
    P_twoPeak=1.0-P_onePeak
    exp_motifOnePeak=round(P_motif*P_onePeak*float(total),0)
    exp_motifTwoPeak=round(P_motif*P_twoPeak*float(total),0)
    exp_nomotifOnePeak=round(P_nomotif*P_onePeak*float(total),0)
    exp_nomotifTwoPeak=round(P_nomotif*P_twoPeak*float(total),0)
    if(P>=0.05): continue
    print(motif,P,sep="\t")
    print("\tmotif\t\tnomotif")
    print("1peak\t"+str(motifOnePeak)+" ("+str(exp_motifOnePeak)+")\t"+
          str(nomotifOnePeak)+" ("+str(exp_nomotifOnePeak)+")")
    print("2peak\t"+str(motifTwoPeaks)+" ("+str(exp_motifTwoPeak)+")\t"+
          str(nomotifTwoPeaks)+" ("+str(exp_nomotifTwoPeak)+")")
    print()



