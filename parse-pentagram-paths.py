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
from SummaryStats import SummaryStats
from Rex import Rex
rex=Rex()

elementsContainingState={}
totalNucsInState={}
transitions={}
twoAndThree=0
totalElements=0
starts={}
ends={}
lengths={}

def processTransitions():
    fromCounts={}
    keys=transitions.keys()
    for trans in keys:
        if(not rex.find("(\d+)->(\d+)",trans)):
            exit("can't parse transition: "+trans)
        fromState=int(rex[1])
        toState=int(rex[2])
        fromCounts[fromState]=fromCounts.get(fromState,0)+transitions[trans]
    for trans in keys:
        rex.find("(\d+)->(\d+)",trans)
        fromState=int(rex[1])
        toState=int(rex[2])
        count=transitions[trans]
        P=float(count)/float(fromCounts[fromState])
        print(str(fromState)+"->"+str(toState)+" = "+str(P))

def parsePath(path):
    parsed=[]
    fields=path.split("|")
    for field in fields:
        if(not rex.find("(\d+):(\d+)-(\d+)",field)):
            exit("error parsing path: "+field)
        state=int(rex[1])
        begin=int(rex[2])
        end=int(rex[3])
        parsed.append([state,begin,end])
    return parsed

def updateStats(parsed):
    global twoAndThree
    containsState=set()
    L=len(parsed)
    firstState=parsed[0][0]
    lastState=parsed[L-1][0]
    starts[firstState]=starts.get(firstState,0)+1
    ends[lastState]=ends.get(lastState,0)+1
    for i in range(L):
        elem=parsed[i]
        (state,begin,end)=elem
        if(lengths.get(state,None) is None): lengths[state]=[]
        lengths[state].append(end-begin)
        containsState.add(state)
        totalNucsInState[state]=totalNucsInState.get(state,0)+end-begin
        if(i+1<L):
            nextState=parsed[i+1][0]
            key=str(state)+"->"+str(nextState)
            transitions[key]=transitions.get(key,0)+1
    for state in containsState:
        elementsContainingState[state]=elementsContainingState.get(state,0)+1
    if(2 in containsState and 3 in containsState): twoAndThree+=1

#=========================================================================
# main()
#=========================================================================
if(len(sys.argv)!=3):
    exit(ProgramName.get()+" <genomewide-pentagram.txt> <LLR-threshold>\n")
(infile,threshold)=sys.argv[1:]
threshold=float(threshold)

IN=open(infile,"rt")
for line in IN:
    fields=line.rstrip().split()
    if(len(fields)!=8): continue
    (fastb,b,e,taskID,LLR,posterior,path,unknown)=fields
    if(not rex.find("t05",fastb)): continue
    if(float(LLR)<threshold): continue
    parsed=parsePath(path)
    updateStats(parsed)
    totalElements+=1
IN.close()

print("\nnumber of elements containing states 2 and 3: "+str(twoAndThree))
print("\ntotal number of elements:"+str(totalElements))
print("\nelements containing each state:")
for state in elementsContainingState.keys():
    count=elementsContainingState[state]
    print("\t"+str(state)+"\t"+str(count))
print("\ntotal nucleotides in each state:")
for state in totalNucsInState:
    count=totalNucsInState[state]
    print("\t"+str(state)+"\t"+str(count))
print("\ntransitions:")
processTransitions()
print("start state:")
for state in starts:
    print("\t"+str(state)+" = "+str(float(starts[state])/float(totalElements)))
print("end state:")
for state in ends:
    print("\t"+str(state)+" = "+str(float(ends[state])/float(totalElements)))
print("state durations:")
for state in lengths.keys():
    [mean,SD,min,max]=SummaryStats.roundedSummaryStats(lengths[state])
    print("\t"+str(state)+"\t"+str(mean)+" +- "+str(SD))
