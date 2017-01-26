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
import os
import sys
import math
from BedReader import BedReader
from Pipe import Pipe
import SumLogProbs
from Interval import Interval
import ProgramName
import TempFilename
from Rex import Rex
rex=Rex()
MUMMIE=os.environ["MUMMIE"]

def getLL(fastb,hmm):
    cmd=MUMMIE+"/get-likelihood "+hmm+" "+fastb
    pipe=Pipe(cmd)
    line=pipe.readline()
    return float(line)

def getPosterior(fastb):
    fgLL=getLL(fastb,fgHMM)
    bgLL=getLL(fastb,bgHMM)
    jointBG=bgLL+bgPrior
    jointFG=fgLL+fgPrior
    denom=SumLogProbs.sumLogProbs2(jointBG,jointFG)
    posterior=math.exp(jointFG-denom)
    return posterior

def getLLR(fastb):
    fgLL=getLL(fastb,fgHMM)
    bgLL=getLL(fastb,bgHMM)
    return fgLL-bgLL

def getPath(hmm,fastb):
    array=[]
    cmd=MUMMIE+"/parse "+hmm+" "+fastb
    pipe=Pipe(cmd)
    pipe.readline() # header
    while(True):
        line=pipe.readline()
        if(not line): break
        if(not rex.find("\S+",line)): continue
        state=int(line)
        array.append(state)    
    return array

def getSections(path):
    sections=[]
    L=len(path)
    elem=[]
    for i in range(L):
        state=path[i]
        eLen=len(elem)
        if(eLen==0 or elem[eLen-1]==state): elem.append(state)
        else:
            sections.append(elem)
            elem=[state]
    if(len(elem)>0): sections.append(elem)
    return sections

def getForeground(path):
    sections=getSections(path)
    begin=None
    end=None
    pos=0
    for section in sections:
        L=len(section)
        type=section[0]
        if(type>1 and type<5 or type>6 and type<10 or type>11 and type<15):
            if(begin is None): begin=pos
        elif(type==5 or type==10 or type==15): end=pos
        pos+=L
    if(end is None): end=pos
    return (begin,end)

def process(dir,posHMM,negHMM,label):
  files=os.listdir()
  n=len(files)
  for i in range(n):
    file=files[i]
    file=file.rstrip()
    if(not rex.find(".fastb",file)): continue
    numer=getLL(file,posHMM)
    denom=getLL(file,negHMM)
    ratio=numer-denom
    ROC.write(str(ratio)+"\t"+str(label)+"\n")

def getParse(path):
    sections=getSections(path)
    if(len(sections)!=5):
        print(path)
        print(sys.argv)
        exit("cannot find five states in parse")
    parse=""
    pos=0
    for section in sections:
        L=len(section)
        type=section[0]
        if(type>1 and type<5 or type>6 and type<10 or type>11 and type<15):
            if(len(parse)>0): parse+=":"
            begin=pos
            end=pos+L
            parse+=str(begin)+"-"+str(end)
        pos+=L
    return parse

#=========================================================================
# main()
#=========================================================================
if(len(sys.argv)!=6):
    exit(ProgramName.get()+
         " <list.txt> <fastb-dir> <pos.hmm> <neg.hmm> <taskID>\n")
(fileList,inDir,fgHMM,bgHMM,taskID)=sys.argv[1:]
files=[]
with open(fileList,"rt") as IN:
    for line in IN:
        files.append(line.rstrip())
nextID=1
for file in files:
    fullPath=inDir+"/"+file
    path=getPath(fgHMM,fullPath)
    whichPath=None;
    if(path[0]==1): whichPath="top"
    elif(path[0]==6): whichPath="middle"
    elif(path[0]==11): whichPath="bottom"
    else: exit("can't identify path")
    (begin,end)=getForeground(path)
    if(begin is None or end is None):
        exit("bad path: "+str(path)+"\n"+fullPath+"\n")
    parse=getParse(path)
    L=end-begin
    llr=getLLR(fullPath)
    id="task"+taskID+"_elem"+str(nextID)
    print(file+"\t"+str(begin)+"\t"+str(end)+"\t"+id+"\t"+str(llr)+"\t"+parse
          +"\t"+whichPath,flush=True)
    nextID+=1

#########################################################################
