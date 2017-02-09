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
import os
import sys
import ProgramName
from BedReader import BedReader
from Rex import Rex
rex=Rex()

#=========================================================================
# main()
#=========================================================================
if(len(sys.argv)!=3):
    exit(ProgramName.get()+" <macs.bed> <t3-1000>\n")
(macsBed,dir)=sys.argv[1:]

keep=set()
files=os.listdir(dir)
for file in files:
    id=file
    if(rex.find("(\S+).t\d+.fastb",file)):
        id=rex[1]
    keep.add(id)

reader=BedReader(macsBed)
while(True):
    record=reader.nextRecord()
    if(record is None): break
    if(not record.chr in keep): continue
    print(record.toString())
reader.close()

