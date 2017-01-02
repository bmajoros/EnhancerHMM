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
from Pipe import Pipe
from Rex import Rex
rex=Rex()

def fisher(a,b,c,d):
    cmd="/home/bmajoros/src/scripts/fisher-exact-test.R "+str(a)+" "+\
        str(b)+" "+str(c)+" "+str(d)
    #pipe=Pipe("/home/bmajoros/src/scripts/fisher-exact-test.R "+str(a)+" "+
    #          str(b)+" "+str(c)+" "+str(d))
    #line=pipe.readline()
    (ret,line)=commands.getstatusoutput(cmd)
    return float(line.rstrip())

if(len(sys.argv)!=4):
    exit(ProgramName.get()+" <hits.txt> <cluster-membership.txt> <cluster-sizes.txt>\n")
(hitsFile,membershipFile,clusterSizeFile)=sys.argv[1:]

clusterSizes={}
IN=open(clusterSizeFile,"rt")
i=0
for line in IN:
    clusterSizes[i]=int(line)
    i+=1
IN.close()

peakToCluster={}
IN=open(membershipFile,"rt")
for line in IN:
    if(not rex.find("\"(\S+)\"\t(\d+)",line)): continue
    peak=rex[1]
    cluster=int(rex[2])
    peakToCluster[peak]=cluster
IN.close()

counts=[0]*10
IN=open(hitsFile,"rt")
for line in IN:
    line=line.rstrip()
    if(not peakToCluster.get(line,None)): continue
    cluster=peakToCluster[line]-1
    counts[cluster]+=1
IN.close()

Pvalues=[]
for i in range(10):
    a=counts[i]
    b=clusterSizes[i]-a
    c=0
    d=0
    for j in range(10):
        if(j==i): continue
        c+=counts[j]
        d+=clusterSizes[j]-counts[j]
    #if(float(a)/float(a+b)<=float(c)/float(c+d)): continue
    P=fisher(a,b,c,d)
    #if(P>=0.05): continue
    proportion=round(float(a)/float(clusterSizes[i]),2)
    expected=round(float(c)/float(c+d),2)
    enrichment=proportion-expected
    if(enrichment<0):
        enrichment=0
        P=1
    msg="cluster "+str(i+1)+"\tP="+str(P)+"\tproportion="+str(proportion)+\
        " vs "+str(expected)+"\t"+str(enrichment)
    Pvalues.append([P,msg])
    #+"\t"+str(float(a)/float(a+b))+" "+str(float(c)/float(c+d)))
    #Pvalues.append([i+1,P,proportion])

#Pvalues.sort(key=lambda x: x[0])
for pair in Pvalues:
    (P,msg)=pair
    print(msg)
