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

cov1=pd.read_csv("cov1.txt",sep=" ",header=None)
cov2=pd.read_csv("cov2.txt",sep=" ",header=None)
cov3=pd.read_csv("cov3.txt",sep=" ",header=None)

max1=max(cov1.max())
max2=max(cov2.max())
max3=max(cov3.max())
m=max(max1,max2,max3)

plt.imshow(cov1, cmap=cm.Reds, interpolation='nearest', vmin=0, vmax=m)
plt.savefig("cov1.pdf")

plt.imshow(cov2, cmap=cm.Reds, interpolation='nearest', vmin=0, vmax=m)
plt.savefig("cov2.pdf")

print(cov3)

plt.imshow(cov3, cmap=cm.Reds, interpolation='nearest', vmin=0, vmax=m)
plt.colorbar()  
plt.savefig("cov3.pdf")


