#!/bin/env python
#=========================================================================
# This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
# License (GPL) version 3, as described at www.opensource.org.
# Copyright (C)2016 William H. Majoros (martiandna@gmail.com).
#=========================================================================
from __future__ import (absolute_import, division, print_function, 
   unicode_literals, generators, nested_scopes, with_statement)
from builtins import (bytes, dict, int, list, object, range, str, ascii,
   chr, hex, input, next, oct, open, pow, round, super, filter, map, zip)
from SlurmWriter import SlurmWriter
import os
import re
from Rex import Rex
rex=Rex()

BASE="/home/bmajoros/GGR/p300"
SLURM_DIR=BASE+"/slurms/max-slurms"
PROGRAM=BASE+"/src/p300-get-p300-max.py"
FASTB=BASE+"/no-dna"
PARTITIONS=BASE+"/partitions"

writer=SlurmWriter()
lists=os.listdir(PARTITIONS)
for list in lists:
  list=list.rstrip()
  if(not rex.find("(\d+)\.txt$",list)): continue
  ID=rex[1]
  infile=PARTITIONS+"/"+list
  outfile=BASE+"/max/max-"+ID+".txt"
  writer.addCommand(PROGRAM+" "+PARTITIONS+"/"+list+" > "+outfile)
writer.mem(5000)
writer.nice(500)
writer.setQueue("new,all")
writer.writeArrayScript(SLURM_DIR,"MAX",SLURM_DIR,500)



