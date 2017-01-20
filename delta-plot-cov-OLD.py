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

if(len(sys.argv)!=3):
    exit(ProgramName.get()+" <in.hmm> <out.means>\n")
(hmmFile,meansfile)=sys.argv[1:]

colormap=cm.bwr
#colormap=cm.Reds

cov1=pd.read_csv("cov1.txt",sep="\t",header=None)
cov2=pd.read_csv("cov2.txt",sep="\t",header=None)
cov3=pd.read_csv("cov3.txt",sep="\t",header=None)

max1=max(cov1.max())
max2=max(cov2.max())
max3=max(cov3.max())
maxAll=max(max1,max2,max3)

#min1=min(cov1.min())
#min2=min(cov2.min())
#min3=min(cov3.min())
#minAll=min(min1,min2,min3)
minAll=-maxAll

plt.imshow(cov1,cmap=colormap,interpolation='nearest',vmin=minAll,vmax=maxAll)
plt.savefig("cov1.pdf")

plt.imshow(cov2,cmap=colormap,interpolation='nearest',vmin=minAll,vmax=maxAll)
plt.savefig("cov2.pdf")

#print(cov3)

plt.imshow(cov3,cmap=colormap,interpolation='nearest',vmin=minAll,vmax=maxAll)
plt.colorbar()  
plt.savefig("cov3.pdf")


