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
from Rex import Rex
rex=Rex()

MOTIFS=("AP1","CEBP","FOX","CTCF","KLF","GR")
BASE="/home/bmajoros/GGR/delta"
SLURM_DIR=BASE+"/slurms/motif-slurms"
PROGRAM=BASE+"/src/motif-scan2.pl"
FASTB=BASE+"/raw"

writer=SlurmWriter()
lists=os.listdir(BASE+"/subsets");
for list in lists:
  list=list.rstrip()
  if(not rex.find("files(\d+)\.txt$",list)): raise Exception(list)
  ID=rex[1]
  for motif in MOTIFS:
    writer.addCommand("cd "+BASE+"\n"+PROGRAM+" "+FASTB+" ../p300/"+
                      motif+" "+motif+
                      " subsets/"+list+" meme/"+ID+"-"+motif+".txt")
writer.mem(5000)
writer.setQueue("new,all")
writer.writeArrayScript(SLURM_DIR,"MEME",1000)



