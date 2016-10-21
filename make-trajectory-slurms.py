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

BASE="/home/bmajoros/GGR/p300"
SLURM_DIR=BASE+"/slurms/trajectory-slurms"
PROGRAM=BASE+"/src/p300-trajectory-posteriors.pl"
FASTB=BASE+"/no-dna"
POS_HMM=BASE+"/model/bootstrap-pos5.hmm"
#NEG_HMM=BASE+"/model/bootstrap-neg5.hmm"
NEG_HMM=BASE+"/model/bootstrap-neg1.hmm"

writer=SlurmWriter()
lists=os.listdir(BASE+"/subsets");
for list in lists:
  list=list.rstrip()
  m=re.search("files(\d+)\.txt$",list)
  if(not m): continue
  ID=m.group(1)
  infile=BASE+"/peak-subsets/"+list
  outfile=BASE+"/trajectories/"+ID+".txt"
  writer.addCommand(PROGRAM+" "+POS_HMM+" "+NEG_HMM+" "+infile+" "+FASTB
                    +" > "+outfile)
writer.mem(5000)
writer.setQueue("new,all")
writer.writeArrayScript(SLURM_DIR,"TRAJECTORIES",SLURM_DIR,500)



