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
import sys
import ProgramName
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import cm

if(len(sys.argv)!=4):
    exit(ProgramName.get()+" <in.cov> <max> <out.pdf>\n")
(infile,maxValue,outfile)=sys.argv[1:]
maxValue=float(maxValue)
minValue=-maxValue

colormap=cm.bwr
cov=pd.read_csv(infile,sep="\t",header=None)
plt.imshow(cov,cmap=colormap,interpolation='nearest',vmin=minValue,
           vmax=maxValue)
plt.colorbar()  
plt.savefig(outfile)


