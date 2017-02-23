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
import ProgramName
from Fastb import Fastb
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
mpl.rcParams['figure.figsize'] = (6,2)
import matplotlib.pyplot as plt
from Rex import Rex
rex=Rex()

#=========================================================================
# main()
#=========================================================================
if(len(sys.argv)!=3):
    exit(ProgramName.get()+" <in.fastb> <out.pdf>\n")
(infile,outfile)=sys.argv[1:]

fastb=Fastb(infile)
for i in range(fastb.numTracks()):
    track=fastb.getIthTrack(i)
    L=track.getLength()
    X=[]
    for j in range(L): X.append(float(j))
    c="black"
    if(rex.find("DNase",track.id)): c="red"
    if(rex.find("P300",track.id)): c="orange"
    if(rex.find("H3K27ac",track.id)): c="blue"
    if(rex.find("H3K4me1",track.id)): c="cyan"
    if(rex.find("H3K4me2",track.id)): c="black"
    plt.plot(X,track.data,color=c)
plt.ylim(-2,2)
plt.savefig(outfile)


