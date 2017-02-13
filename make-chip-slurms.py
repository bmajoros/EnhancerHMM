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

BASE="/home/bmajoros/GGR/delta"
DELTA=BASE
SLURM_DIR=BASE+"/slurms/chip-slurms"
TIMEPOINTS=("t05","t1","t2","t3","t4","t5","t6","t7","t8","t10","t12")

writer=SlurmWriter()
for time in TIMEPOINTS:
  writer.addCommand("cd "+DELTA+" ; src/analyze-raw-p300.py genomewide-loopy.txt    1000 250 "+time+" > chip/loopy-p300-"+time+".txt")
  writer.addCommand("cd "+DELTA+" ; src/analyze-raw-GR.py genomewide-loopy.txt      1000 250 "+time+" > chip/loopy-gr-"+time+".txt")
  writer.addCommand("cd "+DELTA+" ; src/analyze-raw-cebp.py genomewide-loopy.txt    1000 250 "+time+" > chip/loopy-cebp-"+time+".txt")
  writer.addCommand("cd "+DELTA+" ; src/analyze-raw-ap1.py genomewide-loopy.txt     1000 250 "+time+" > chip/loopy-ap1-"+time+".txt")
  writer.addCommand("cd "+DELTA+" ; src/analyze-raw-dnase.py genomewide-loopy.txt   1000 250 "+time+" > chip/loopy-dnase-"+time+".txt")
  writer.addCommand("cd "+DELTA+" ; src/analyze-raw-h3k4me1.py genomewide-loopy.txt 1000 250 "+time+" > chip/loopy-h3k4me1-"+time+".txt")
  writer.addCommand("cd "+DELTA+" ; src/analyze-raw-h3k4me2.py genomewide-loopy.txt 1000 250 "+time+" > chip/loopy-h3k4me2-"+time+".txt")
  writer.addCommand("cd "+DELTA+" ; src/analyze-raw-h3k27ac.py genomewide-loopy.txt 1000 250 "+time+" > chip/loopy-h3k27ac-"+time+".txt")
writer.mem(5000)
writer.nice(500)
writer.setQueue("new,all")
writer.writeArrayScript(SLURM_DIR,"CHIP",SLURM_DIR,500)



