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

#CLASSES="0,2,3,2,3,2,1,2,3,2,1"
#CLASSES="0,1,2,3,2,1,1,2,3,2,1,2,3,2,3,2"
CLASSES="0,1,2,3,2,1"
MUMMIE=os.environ["MUMMIE"]
BASE="/home/bmajoros/GGR/delta"
SLURM_DIR=BASE+"/slurms/annotate-slurms"
PROGRAM=BASE+"/src/annotate-partition.pl"
FASTB=BASE+"/delta-fastb"
POS_HMM=BASE+"/hmm/loopy-1path-trained2.hmm"
PARTITIONS=BASE+"/fastb-partitions"
OUTDIR=BASE+"/annotated-loopy"

writer=SlurmWriter()
lists=os.listdir(PARTITIONS)
for list in lists:
  list=list.rstrip()
  if(not rex.find("(\d+)\.txt$",list)): continue
  ID=rex[1]
  infile=PARTITIONS+"/"+list
  writer.addCommand(PROGRAM+" "+FASTB+" "+POS_HMM+" "+CLASSES
                    +" "+infile+" "+OUTDIR)
writer.mem(5000)
writer.nice(500)
writer.setQueue("new,all")
writer.writeArrayScript(SLURM_DIR,"ANNOTATE",SLURM_DIR,500)



