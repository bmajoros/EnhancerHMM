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

if(len(sys.argv)!=3):
    exit(sys.argv[0]+" <in.sam> <out.sam>")
(infile,outfile)=sys.argv[1:]

with open(outfile,"wt") as OUT:
    with open(infile,"rt") as IN:
        for line in IN:
            if(len(line)<1): continue
            if(line[0]=="@"):
                OUT.write(line)
                continue
            fields=line.split()
            if(fields[5]!="25M" and fields[5]!="26M"): continue
            OUT.write(line)



