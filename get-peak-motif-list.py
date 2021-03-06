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
import os
from Rex import Rex
rex=Rex()

def process(dir,motif):
    files=os.listdir(dir)
    for file in files:
        if(rex.find("(\S+)\.standardized",file)):
            peak=rex[1]
            print(peak,motif,sep="\t")

BASE="/home/bmajoros/GGR/delta"
process("CTCF-fastb","CTCF")
process("FOX-fastb","FOX")
process("CEBP-fastb","CEBP")
process("AP1-fastb","AP1")
process("KLF-fastb","KLF")
process("GR-fastb","GR")



