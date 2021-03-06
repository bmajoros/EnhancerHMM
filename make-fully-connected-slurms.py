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
import sys
import os
import random
from SlurmWriter import SlurmWriter

REPLICATES=100
BASE="/home/bmajoros/GGR/delta/fastb"
TRAINING_SETS_POS=BASE+"/positives"
MUMMIE=os.environ["MUMMIE"]
BAUMWELCH=MUMMIE+"/baum-welch"
HMM_DIR="/home/bmajoros/GGR/delta/hmm"
TEMPLATE_HMM=HMM_DIR+"/fully-connected.hmm"
NEG_TEMPLATE_HMM=HMM_DIR+"/trained-neg-nomotif.hmm"
SLURM_DIR="/home/bmajoros/GGR/delta/slurms/pentagram-slurms"

#=========================================================================
# main()
#=========================================================================

writer=SlurmWriter()
for rep in range(REPLICATES):
    outfile=HMM_DIR+"/pentagram"+"-rep"+str(rep)+".hmm"
    cmd="cd "+HMM_DIR+" ; "\
        +BAUMWELCH\
        +" -s "+str(random.randint(0,1000000000))\
        +" -S -W -t tie.txt -c 32 -L 0.01 -N 10000 "\
        +TEMPLATE_HMM+" tgf.tgf "+TRAINING_SETS_POS+" 1000 "+outfile
    writer.addCommand(cmd)
writer.addCommand(cmd)
writer.mem(5000)
writer.nice(500)
writer.threads(32)
writer.setQueue("new,all")
writer.writeArrayScript(SLURM_DIR,"PENTAGRAM",SLURM_DIR,500)


