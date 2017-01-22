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

NUM_JOBS=100
BASE="/home/bmajoros/GGR/delta"
SLURM_DIR=BASE+"/slurms/train-slurms"
PROGRAM="/home/bmajoros/src/MUMMIE/baum-welch"
#trainDir="../train-pos-nomotif"
trainDir="../train-neg-nomotif"

writer=SlurmWriter()
for i in range(NUM_JOBS):
  #outfile="twopaths-trained"+str(i)+".hmm"
  outfile="neg-2paths"+str(i)+".hmm"
  writer.addCommand("cd /home/bmajoros/GGR/delta/hmm ; "+
                    PROGRAM+" -W -t tie.txt -c 32 -L 0.01 -N 1000 twopaths.hmm"
                    +" tgf.tgf "+trainDir+" 1000 "+outfile)
writer.mem(5000)
writer.nice(500)
writer.threads(32)
writer.setQueue("new,all")
writer.writeArrayScript(SLURM_DIR,"EM",SLURM_DIR,500)


