import os
import sys
import random
import time
import math
import copy
import re
import itertools
from collections import defaultdict


def closure(a):
    while True:
        changed=False
        for i,x in enumerate(a):
            for j,y in enumerate(a):
                if x==y:
                    continue
                tmp = list(set(x+y))
                if len(tmp)<(len(x)+len(y)):
                    a.append(tmp)
                    a.remove(x)
                    a.remove(y)
                    changed=True
                    break
            if changed:
                break
        if not changed:
            break
    return a
    
def generateFeatures(mlnfile,cbound,jbound):
    basefile = os.path.basename(mlnfile)
    tmpspecfile = "tmp-"+basefile[:basefile.find(".")]+".spec"
    tmpgmfile = "tmp-"+basefile[:basefile.find(".")]+".gm"
    ifile = open(mlnfile)
    first = ifile.readline()
    predparts = first.split()
    preds = []
    domsizes = []
    cnt = 0
    for p in predparts:
        parts1 = p.strip().split(":")
        if len(parts1)>3:
            sys.exit("Currently only supports unary and binary predicates")
        preds.append(parts1[0])
        domsizes.append(parts1[1:])
    alljoins = []
    mlines = ifile.readlines()
    ifile.close()
    tmpmln = []
    for line in mlines:
        if len(line) < 2:
            continue
        tmpmln.append(line)
        parts = line.split(":")
        atoms = parts[1].split(" v ")
        #print(atoms)
        variables = []
        prednames = []
        tmpset = []
        for a in atoms:
            ind = a.find("(")
            predname = a[:ind]
            if predname[0]=="!":
                predname = predname[1:]
            prednames.append(predname)
            rest = a[ind+1:]
            rest = rest[:rest.find(")")]
            varnames = rest.split(",")
            variables.append(varnames)
            for v in varnames:
                if v not in tmpset:
                    tmpset.append(v)
        for t in tmpset:
            joins = []
            for j,p in enumerate(prednames):
                for ix,k in enumerate(variables[j]):
                    if t==k and (p+"-id"+str(ix)) not in joins:
                        joins.append(p+"-id"+str(ix))
            alljoins.append(joins)
    #print(alljoins)
    closure(alljoins)
    alljoins1 = []
    for a in alljoins:
        if a not in alljoins1:
            alljoins1.append(a)
    alljoins = alljoins1
    #print("****")
    #print(alljoins)       
    ofile = open(tmpspecfile,'w')
    ofile.write(first)
    ofile.write(str(len(alljoins))+"\n")
    for jnlist in alljoins:
        tmp = jnlist[0].split("-")
        pidx = preds.index(tmp[0])
        idx = int(tmp[1][2:])
        newdom = int(float(cbound)*float(domsizes[pidx][idx]))
        if newdom==0:
            newdom=1
        ofile.write(str(newdom)+":"+",".join(jnlist)+"\n")
 
    for i,p in enumerate(preds):
        ofile.write(p+":")
        tmpstr = ""
        for j in range(0,len(domsizes[i]),1):
            tmp = p+"-id"+str(j)
            for k,a in enumerate(alljoins):
                if tmp in a:
                    tmpstr= tmpstr + str(k)+",0:"
                    break
        ofile.write(tmpstr[:len(tmpstr)-1]+" ")
    ofile.write("\n"+str(len(tmpmln))+"\n")
    for t in tmpmln:
        ofile.write(t)
    ofile.close()

generateFeatures(sys.argv[1],sys.argv[2],1)
