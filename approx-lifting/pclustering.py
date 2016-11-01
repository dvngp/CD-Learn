import os
import sys
import math
import itertools
import numpy as np
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
import pandas as pd

def calc_dist(p1,p2):
    res = 0.0
    for i in range(0,len(p1),1):
        res = res + ((p2[0] - p1[0]) ** 2)
    return math.sqrt(res ** 2)
    
def intersect(a,b):
    return list(set(a) & set(b))

def doclustering(mlnfile,evidfile,cratio,outmfile,outefile,outclog):
    #exepath = "./Release/cdlearn"
    basefile = os.path.basename(mlnfile)
    os.system('python genfeatures.py '+mlnfile+" "+cratio+"\n")
    tmpspecfile  = "tmp-"+basefile[:basefile.rfind(".")]+".spec"

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

    ifile = open(evidfile)
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
            predtables[tidx][i][idx] = predtables[tidx][i][idx]+5
    ifile.close()
    #print(predtables)


    ifile = open(evidfile)
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
    #sys.exit(0)
    allclusters = []
    numclusters = []
    for r in reducedsizes:
        numclusters.append(r)
    for i,t in enumerate(tables):
        idx = prednames.index(t[0])
        sz = len(predtables[idx][tableids[i][0]])
        features = []
        for j in range(0,sz,1):
            v_features = []
            for k in range(0,len(t),1):
                ix = prednames.index(t[k])
                v_features.append(predtables[ix][tableids[i][k]][j])
                ix1 = pvec.index(t[k])
                #print(str(float(predtables[ix][tableids[i][k]][j])/sizevec[ix1][tableids[i][k]]))
            
            features.append(v_features)
        clustervec = []
        #print("numclusters="+str(numclusters[i]))
        for j in range(0,numclusters[i],1):
            tmp = []
            clustervec.append(tmp)
        #if os.name=="nt":
        #    os.system("java -cp .;./weka.jar JCLift "+"d-"+str(i)+".csv "+sys.argv[5]+" "+str(numclusters[i])+" 102923810 > d-"+str(i)+".dat")
        #else:
        #    os.system("java -cp .:./weka.jar JCLift "+"d-"+str(i)+".csv "+sys.argv[5]+" "+str(numclusters[i])+" 102923810 > d-"+str(i)+".dat")
        
        #print(Z[i])
        
        print("Clustering using KMeans")
        KM = KMeans(n_clusters=numclusters[i])
        KM.fit(features)
            
        '''
        ifile = open("d-"+str(i)+".dat")
        lines = ifile.readlines()
        ifile.close()
        intvals = []
        for j,l in enumerate(lines):
            clustervec[int(l.strip())].append(j)
        '''
        for j,l in enumerate(KM.labels_):
            clustervec[l].append(j)
        allclusters.append(clustervec)
    #print(allclusters)

    print("Processing EVidence...")        
    ofile=open(outefile,'w')
    numcomponents = 0
    for i,p in enumerate(pvec):
        ix = prednames.index(p)
        #ofile1 = open(p+".map",'w')
        if len(pindvec[i])==1:
            for j in range(0,len(allclusters[pindvec[i][0]]),1):
                if len(allclusters[pindvec[i][0]][j])==0:
                    ofile.write("0 "+p+"("+str(j)+")\n")
                else:
                    v = len(intersect(sevidtable[ix],allclusters[pindvec[i][0]][j]))/float(len(allclusters[pindvec[i][0]][j]))
                    numcomponents = numcomponents + len(allclusters[pindvec[i][0]][j])
                    ofile.write(str(v)+" "+p+"("+str(j)+")\n")
                #ofile1.write(str(len(allclusters[pindvec[i][0]][j]))+"\n")
        elif len(pindvec[i])==2:
            for j in range(0,len(allclusters[pindvec[i][0]]),1):
                for k in range(0,len(allclusters[pindvec[i][1]]),1):  
                    count = 0 
                    total = 0
                    tmplist = []
                    tmplist.append(allclusters[pindvec[i][0]][j]) 
                    tmplist.append(allclusters[pindvec[i][1]][k]) 
                    for px in itertools.product(*tmplist):
                        total = total+1
                        for m in range(0,len(evidtable[ix][px[0]]),1):
                            if evidtable[ix][px[0]][m]==px[1]:
                                count = count + 1
                                break
                    #print(str(total)+" "+str(count))
                    #print(str(j)+" "+str(k))
                    if total==0:
                        ofile.write("0 "+p+"("+str(j)+","+str(k)+")\n")
                        numcomponents = numcomponents + 1
                        
                    else:
                        ofile.write(str(count/float(total))+" "+p+"("+str(j)+","+str(k)+")\n")
                        numcomponents = numcomponents + 1
                    #ofile1.write(str(total)+"\n")
        #ofile1.close()
    ofile.close()

    '''
    #Write compressed MLN in Alchemy format
    ofile=open("cmln.txt",'w')
    doms = []
    preddefs = []
    for i,t in enumerate(tables):
        doms.append("dom"+str(i)+"={0, ..., "+str(numclusters[i])+"}\n")
    for i,p in enumerate(pvec):
        ix = prednames.index(p)
        s = prednames[ix]+"("
        for j,l in enumerate(pindvec[ix]):
            s = s+"dom"+str(l)
            if j!=len(pindvec[ix])-1:
                s = s + ","
        s = s + ")\n"
        preddefs.append(s)
    for d in doms:
        ofile.write(d)
    for p in preddefs:
        ofile.write(p)
    for o in otherlines:
        p = o.split(":")
        ofile.write(p[0]+" "+p[1])
    ofile.close()
    '''

    doms = []
    for i,t in enumerate(tables):
        doms.append(numclusters[i]) 
    ofile=open(outmfile,'w')
    for i,p in enumerate(pvec):
        ix = prednames.index(p)
        ofile.write(prednames[ix]+":")
        for j,l in enumerate(pindvec[ix]):
            ofile.write(str(doms[l]))
            if j!=len(pindvec[ix])-1:
                ofile.write(":")
        ofile.write(" ")
    ofile.write("\n")
    for o in otherlines:
        ofile.write(o)
    ofile.close()

    #Log the clusters
    ofile = open(outclog,'w')
    for i in range(0,len(allclusters),1):
        for j in range(0,len(allclusters[i]),1):
            ofile.write(str(allclusters[i][j])+" ")
        ofile.write("\n")
    ofile.close()


    '''
    print("Performing Inference...")
    #os.system("./cdlearn cmln.txt cevid.txt cquery.txt 10000 1 MAR res.dat")
    #os.system(exepath+" cmln.txt cevid.txt " + sys.argv[3]+ " "+str(numcomponents*10) + " 1 MAR res.dat "+" "+sys.argv[7])
    os.system(exepath+" cmln.txt cevid.txt " + sys.argv[3]+ " 30 1 MAR res.dat dum.dat dum.dat")
    print("Collecting Results...")
    #project back the results to original dimensions
    ifile=open("res.dat")
    lines = ifile.readlines()
    ifile.close()
    ofile = open("margs.dat",'w')
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
    '''