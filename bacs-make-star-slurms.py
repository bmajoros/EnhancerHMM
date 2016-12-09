#!/usr/bin/env python
#=========================================================================
# This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
# License (GPL) version 3, as described at www.opensource.org.
# Copyright (C)2016 William H. Majoros (martiandna@gmail.com).
#=========================================================================
from __future__ import (absolute_import, division, print_function, 
   unicode_literals, generators, nested_scopes, with_statement)
from builtins import (bytes, dict, int, list, object, range, str, ascii,
   chr, hex, input, next, oct, open, pow, round, super, filter, map, zip)
# The above imports should allow this program to run in both Python 2 and
# Python 3.  You might need to update your version of module "future".
from SlurmWriter import SlurmWriter
import glob
from Rex import Rex
rex=Rex()

# Global variables
STAR="/data/reddylab/software/STAR_2.4.2a/STAR-STAR_2.4.2a/bin/Linux_x86_64/STAR";
numThreads=8
maxParallel=500
memory=40000
jobName="STAR"
BASE="/home/bmajoros/GGR/bacs"
#dnaFastqDir="/gpfs/fs0/data/reddylab/Flowcells/Project_Vockley_709_131016"
#rnaFastqDir1="/data/reddylab/Flowcells/Project_Vockley_743_131102"
#rnaFastqDir2="/data/reddylab/Flowcells/Project_Vockley_744_131108"
fastqDir=BASE+"/fastq"
slurmDir=BASE+"/star-slurms"
starIndex=BASE+"/STAR"
outputDir=starIndex+"/output"
slurm=SlurmWriter()

def process(fastqDir,tag):
    files=glob.glob(fastqDir+"/*.fastq")
    for file in files:
        if(not rex.find("(\S+)Read1.fastq",file)): continue
        pair=rex[1]+"Read2.fastq"
        if(not rex.find("([^/]+)[_-]Read1.fastq",file)): raise Exception("")
        filestem=rex[1]+tag
        cmd=STAR+" --genomeLoad LoadAndKeep --genomeDir "+starIndex+" --readFilesIn "+file+" "+pair+" --outFileNamePrefix "+filestem+" --runThreadN "+str(numThreads)
        slurm.addCommand(cmd)

#=========================================================================
#                               main()
#=========================================================================
#process(dnaFastqDir,"")
#process(rnaFastqDir1,"_x1")
#process(rnaFastqDir2,"_x2")
process(fastqDir,"")
slurm.mem(memory)
slurm.threads(numThreads)
slurm.setQueue("new,all")
slurm.writeArrayScript(slurmDir,jobName,starIndex,maxParallel)






