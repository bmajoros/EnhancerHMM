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

NUM_PARTITIONS=5
REPLICATES=100
BASE="/home/bmajoros/GGR/delta/fastb"
CROSSVAL=BASE+"/crossval"
TRAINING_SETS=CROSSVAL+"/training-sets"
TRAINING_SETS_POS=TRAINING_SETS+"/pos"
TRAINING_SETS_NEG=TRAINING_SETS+"/neg"
MUMMIE=os.environ["MUMMIE"]
BAUMWELCH=MUMMIE+"/baum-welch"
HMM_DIR="/home/bmajoros/GGR/delta/hmm"
TEMPLATE_HMM=HMM_DIR+"/fully-connected.hmm"
NEG_TEMPLATE_HMM=HMM_DIR+"/trained-neg-nomotif.hmm"
SLURM_DIR="/home/bmajoros/GGR/delta/slurms/fully-connected-slurms"

#=========================================================================
# main()
#=========================================================================

writer=SlurmWriter()
for i in range(NUM_PARTITIONS):
    bin=i+1
    traindir=TRAINING_SETS_NEG+"/leaveout"+str(bin)
    outfile=HMM_DIR+"/crossval-neg-bin"+str(bin)+".hmm"
    cmd="cd "+HMM_DIR+" ; "+BAUMWELCH+" -s "\
        +str(random.randint(0,1000000000))+" -S -c 32 -L 0.01 -N 10000 "\
        +NEG_TEMPLATE_HMM+" tgf.tgf "+traindir+" 1000 "+outfile
    writer.addCommand(cmd)
    for rep in range(REPLICATES):
        traindir=TRAINING_SETS_POS+"/leaveout"+str(bin)
        outfile=HMM_DIR+"/crossval-pos-bin"+str(bin)+"-rep"+str(rep)+".hmm"
        cmd="cd "+HMM_DIR+" ; "+BAUMWELCH+" -s "\
            +str(random.randint(0,1000000000))\
            +" -S -W -t tie.txt -c 32 -L 0.01 -N 10000 "\
            +TEMPLATE_HMM+" tgf.tgf "+traindir+" 1000 "+outfile
        writer.addCommand(cmd)
writer.mem(5000)
writer.nice(500)
writer.threads(32)
writer.setQueue("new,all")
writer.writeArrayScript(SLURM_DIR,"PENTAGRAM",SLURM_DIR,500)


