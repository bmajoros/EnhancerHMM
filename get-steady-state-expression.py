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

DIR="/data/reddylab/projects/GGR/data/rna_seq/quantified_read_counts/iter0/"
REPS=(1,2,3,4)
NUM_REPS=len(REPS)

byGene={}
for rep in REPS:
    file=DIR+"t00_rep"+str(rep)+".rsem.genes.results"
    with open(file) as IN:
        IN.readline() # header
        for line in IN:
            fields=line.rstrip().split()
            if(len(fields)!=7): continue
            (gene,transcripts,length,effLen,expCount,TPM,FPKM)=fields
            FPKM=float(FPKM)
            byGene[gene]=byGene.get(gene,0.0)+FPKM/float(NUM_REPS)
genes=byGene.keys()
for gene in genes:
    print(gene,byGene[gene],sep="\t")

