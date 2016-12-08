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
rex=Rex()

def setHash(hash,chr,pos,value):
    if(hash.get(chr,None) is None): hash[chr]=[]
    hash[chr][pos]=value

def load(file,hash):
    with open(file,"rt") as fh:
        header=fh.readline()
        for line in fh:
            fields=line.split()
            if(len(fields)!=6): continue
            (label,fc,cpm,lr,p,fdr)=fields
            if(fc<=0 || fdr>0.05): continue
            if(not rex.find("(\S+)_(\d+)$",label)):
                raise Exception("can't parse line: "+line)
            chr=rex[1]; pos=rex[2]
            setHash(hash,chr,pos,fc)

if(len(sys.argv)!=3):
    exit(sys.argv[0]+" <in.eth> <in.dex>")
(ethFile,dexFile)=sya.argv[1:]

hash={}
load(ethFile,hash)
load(dexFile,hash)



