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
    exit(ProgramName.get()+" <clusters.tsv> <cluster-sizes.tsv>\n")
(clusterFile,sizeFile)=sys.argv[1:]

df=pd.read_csv(clusterFile,sep="\t",header=None)
sizes=pd.read_csv(sizeFile,sep="\t",header=None)

#print(sizes)
#print(df)
#df=[x for (y,x) in sorted(zip(sizes,df), key=lambda pair: pair[0])]
#print(df)
#df=pd.DataFrame(df)
#print(df)

#columns=df.columns.values
#columns[0]='a'
#df.columns=columns
#df.sort_values(by='a',inplace=True)

plt.imshow(df, cmap=cm.Reds, interpolation='nearest')
plt.colorbar()  
plt.savefig("heatmap.pdf")

plt.imshow(sizes, cmap=cm.Blues, interpolation='nearest')
plt.colorbar()  
plt.savefig("sizes.pdf")

