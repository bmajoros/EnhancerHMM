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
from scipy.stats.stats import spearmanr
import random
import copy
import numpy as np
from sklearn import linear_model
import statsmodels.api as sm
from scipy import stats

REPS=1000

def load(filename):
    matrix=[]
    with open(filename,"rt") as IN:
        for line in IN:
            fields=line.rstrip().split()
            row=[]
            for field in fields:
                row.append(float(field))
            matrix.append(row)
    return matrix

def compareRows(row1,row2):
    L=len(row1)
    if(len(row2)!=L): exit("rows have unequal lengths")
    (r,p)=spearmanr(row1,row2)
    extreme=0
    for i in range(REPS):
        #permuted=permute(row2)
        permuted=copy.copy(row2)
        for i in range(L): permuted[i]=random.uniform(-1.0,1.0) ###
        (rr,pp)=spearmanr(row2,permuted)
        if(r>0 and rr>=r or r<0 and rr<=r): extreme+=1
    newP=float(extreme)/float(REPS)
    print(r,p,newP,sep="\t")

def permute(row):
    L=len(row)
    newRow=[]
    for a in row: newRow.append(a)
    for i in range(L):
        j=random.randint(i,L-1)
        temp=newRow[i]
        newRow[i]=newRow[j]
        newRow[j]=temp
    return newRow


#=========================================================================
# main()
#=========================================================================
if(len(sys.argv)!=3):
    exit(ProgramName.get()+" <enhancer-scores.txt> <gene-expression-scores.txt\n")
(enhancerFile,geneFile)=sys.argv[1:]

enhancerData=load(enhancerFile)
geneData=load(geneFile)
numRows=len(enhancerData)
if(len(geneData)!=numRows): exit("mismatch in number of rows")
bigRow1=[]; bigRow2=[]
for i in range(numRows):
    #compareRows(enhancerData[i],geneData[i])
    bigRow1.extend(enhancerData[i])
    bigRow2.extend(geneData[i])
compareRows(bigRow1,bigRow2)

bigRow1=np.array(bigRow1)
bigRow2=np.array(bigRow2)
#print(len(bigRow1),len(bigRow2))
bigRow1=np.reshape(bigRow1,(-1,1))
bigRow2=np.reshape(bigRow2,(-1,1))
regr = linear_model.LinearRegression()
regr.fit(bigRow1,bigRow2)

print('Coefficients: \n', regr.coef_)
print("Mean squared error: %.2f"
      % np.mean((regr.predict(bigRow1) - bigRow2) ** 2))
print('Variance score: %.2f' % regr.score(bigRow1, bigRow2))

X=bigRow1
y=bigRow2

X2 = sm.add_constant(X)
est = sm.OLS(y, X2)
est2 = est.fit()
print(est2.summary())

