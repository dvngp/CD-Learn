import os
import sys
import math
import itertools

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


def intersect(a,b):
    return list(set(a) & set(b))

if len(sys.argv)<5:
    sys.exit("Usage: projectresults mlnfile evidfile inresfile clogfile")

basefile = os.path.basename(sys.argv[1])
generateFeatures(sys.argv[1],1,1)
tmpspecfile  = "tmp-"+basefile[:basefile.rfind(".")]+".spec"

#ep1 = float(sys.argv[4])
#ep2 = float(sys.argv[5])

ifile=open(tmpspecfile)
firstline = ifile.readline()
numdoms = int(ifile.readline().strip())
reducedsizes = []
tables = []
tableids = []
for i in range(0,numdoms,1):
    domline = ifile.readline()
    domvals = domline.strip().split(":")
    reducedsizes.append(int(domvals[0]))
    featureids = domvals[1].split(",")
    namevec = []
    valvec = []
    for d in featureids:
        d1 = d.split("-")
        namevec.append(d1[0])
        valvec.append(int(d1[1][2:]))
    tables.append(namevec)
    tableids.append(valvec)
#print(tables)
#print(tableids)        
#print(reducedsizes)

parts = firstline.strip().split()
prednames = []
predtables = []
for p in parts:
    p1 = p.split(":")
    prednames.append(p1[0])
    tmp1 = []
    for j in range(1,len(p1),1):
        sz = int(p1[j])
        tmp = []
        for k in range(0,sz,1):
            tmp.append(0)
        tmp1.append(tmp)
    predtables.append(tmp1)

    
'''
ofile = open("cquery.txt",'w')
for p in prednames:
    ofile.write(p+"\n")
ofile.close()
'''

evidtable = []
sevidtable = []
for i,p in enumerate(parts):
    p1 = p.split(":")
    tmp2 = []
    stmp1 = []
    for j in range(1,len(p1),1):
        sz = int(p1[j])
        for k in range(0,sz,1):
            if len(p1)==2:
                stmp1.append(-1)
            if j==1:
                tmp3 = []
                tmp2.append(tmp3)                
    evidtable.append(tmp2)
    sevidtable.append(stmp1)

parts = ifile.readline().strip().split()
pvec = []
pindvec = []
sizevec = []
for p in parts:
    parts1 = p.split(":")
    predname = parts1[0]
    predindexes = []
    for p1 in parts1[1:]:
        p2 = p1.split(",")
        predindexes.append(int(p2[0]))
    pvec.append(predname)
    pindvec.append(predindexes)
    ix = prednames.index(predname)
    sizes = []
    for i in range(0,len(predtables[ix]),1):
        sz = 1
        for j in range(0,len(predtables[ix]),1):
            if i==j:
                continue
            sz = sz*len(predtables[ix][j])
        sizes.append(sz)
    sizevec.append(sizes)
#print(pvec)
#print(pindvec)
#print(sizevec)

ifile.readline()
otherlines=ifile.readlines()
ifile.close()

ifile = open(sys.argv[2])
for l in ifile:
    parts = l.strip().split("(")
    if parts[0][0]=='!':
        continue
    arglist = parts[1][:parts[1].find(")")].split(",")
    for i,a in enumerate(arglist):
        idx = int(a.strip())
        if idx==-1:
            continue
        tidx = -1
        try:
            tidx = prednames.index(parts[0])
        except ValueError:
            continue
        predtables[tidx][i][idx] = predtables[tidx][i][idx]+1.0
ifile.close()


ifile = open(sys.argv[2])
evidlines = ifile.readlines()
totalevidences = len(evidlines)
for l in evidlines:
    parts = l.strip().split("(")
    ix = -1
    try:
        ix = prednames.index(parts[0])
    except ValueError:
        continue
    arglist = parts[1][:parts[1].find(")")].split(",")
    argsint = []
    for i,a in enumerate(arglist):
        argsint.append(int(a.strip()))
    if argsint[0]==-1:
        continue
    if len(argsint)==1:
        sevidtable[ix][argsint[0]]=argsint[0]   
    for a in argsint[1:]:
        evidtable[ix][argsint[0]].append(a)
#print(evidtable)        
ifile.close()

'''
for i,t in enumerate(tables):
    idx = prednames.index(t[0])
    sz = len(predtables[idx][tableids[i][0]])
    numfeatures = len(t)
    for k in range(0,numfeatures,1):
        ix = prednames.index(t[k])
        v = sum(predtables[ix][tableids[i][k]])
        if v==0:
            continue
        for j in range(0,sz,1):
            predtables[ix][tableids[i][k]][j] = predtables[ix][tableids[i][k]][j]/float(v)
'''
allclusters = []
ifile = open(sys.argv[4])
for l in ifile:
    if len(l) < 2:
        continue
    parts = l.strip().split("]")
    tmp1 = []
    for p1 in parts[:len(parts)-1]:
        iarr = p1.split(",")
        iarr[0] = iarr[0][iarr[0].find("[")+1:]
        #iarr[len(iarr)-1] = iarr[len(iarr)-1][:len(iarr[len(iarr)-1])-1]
        tmp = []
        for ix in iarr:
            if len(ix)==0:
                continue
            tmp.append(int(ix))
        tmp1.append(tmp)
    allclusters.append(tmp1)
ifile.close()
#print(allclusters)
ifile=open(sys.argv[3])
lines = ifile.readlines()
ifile.close()
ofile = open("pr-"+sys.argv[4],'w')
for l in lines:
    if len(l)<2:
        continue
    prt = l.strip().split()
    prob = float(prt[0])
    parts = prt[1].strip().split("(")
    #if parts[0][0]=='!':
    #    continue
    arglist = parts[1][:parts[1].find(")")].split(",")
    tidx = prednames.index(parts[0])
    tidx1 = pvec.index(parts[0])
    idx = int(arglist[0].strip())
    idx1 = pindvec[tidx1][0]
    if len(arglist)==1:
        for v in allclusters[idx1][idx]:
            ofile.write(str(prob)+" "+parts[0]+"("+str(v)+")\n")
    else:
        idx2 = int(arglist[1].strip())
        idx3 = pindvec[tidx1][1]
        tmp = []
        tmp.append(allclusters[idx1][idx])
        
        tmp.append(allclusters[idx3][idx2])

        for v in itertools.product(*tmp):
            st = str(v[0])
            for r in v[1:]:
                st = st +","+ str(r)
            ofile.write(str(prob)+" "+parts[0]+"("+st+")\n")
ofile.close()
