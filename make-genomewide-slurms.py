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
SLURM_DIR=BASE+"/slurms/genomewide-slurms"
PROGRAM=BASE+"/src/run-genome-wide.py"
FASTB=BASE+"/no-dna"
POS_HMM=BASE+"/model/trained-pos5.hmm"
NEG_HMM=BASE+"/model/trained-neg1.hmm"
PARTITIONS=BASE+"/partitions"

writer=SlurmWriter()
lists=os.listdir(PARTITIONS)
for list in lists:
  list=list.rstrip()
  if(not rex.find("(\d+)\.txt$",list)): continue
  ID=rex[1]
  infile=PARTITIONS+"/"+list
  outfile=BASE+"/genomewide-"+ID+".txt"
  writer.addCommand(PROGRAM+" "+PARTITIONS+"/"+list+" "+FASTB+" "
                    +POS_HMM+" "+NEG_HMM+" "+ID+" > "+outfile)
writer.mem(5000)
writer.nice(500)
writer.setQueue("new,all")
writer.writeArrayScript(SLURM_DIR,"GENOME",SLURM_DIR,500)



