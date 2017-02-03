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

BASE="/home/bmajoros/GGR/delta"
SLURM_DIR=BASE+"/slurms/run-genomewide-slurms"
PROGRAM=BASE+"/src/run-genome-wide.py"
FASTB=BASE+"/delta-fastb"
#POS_HMM=BASE+"/hmm/trained-nomotif.hmm"
#POS_HMM=BASE+"/hmm/twopaths-bestLL.hmm"
#POS_HMM=BASE+"/hmm/threepaths-best.hmm"
POS_HMM=BASE+"/hmm/1path-best-full.hmm"
NEG_HMM=BASE+"/hmm/trained-neg-nomotif.hmm"
PARTITIONS=BASE+"/fastb-partitions"

writer=SlurmWriter()
lists=os.listdir(PARTITIONS)
for list in lists:
  list=list.rstrip()
  if(not rex.find("(\d+)\.txt$",list)): continue
  ID=rex[1]
  infile=PARTITIONS+"/"+list
  outfile=BASE+"/genomewide/"+ID+".txt"
  writer.addCommand(PROGRAM+" "+infile+" "+FASTB+" "
                    +POS_HMM+" "+NEG_HMM+" "+ID+" > "+outfile)
writer.mem(5000)
writer.nice(500)
writer.setQueue("new,all")
writer.writeArrayScript(SLURM_DIR,"GENOMEWIDE",SLURM_DIR,500)



