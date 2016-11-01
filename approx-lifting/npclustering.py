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
    
def generateFeatures(mlnfile,cbound,jbound,r):
    basefile = os.path.basename(mlnfile)
    tmpspecfile = str(r)+"-tmp-"+basefile[:basefile.find(".")]+".spec"
    tmpgmfile = str(r)+"-tmp-"+basefile[:basefile.find(".")]+".gm"
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

	
def calc_dist(p1,p2):
    res = 0.0
    for i in range(0,len(p1),1):
        res = res + ((p2[0] - p1[0]) ** 2)
    return math.sqrt(res ** 2)


def intersect(a,b):
    return list(set(a) & set(b))

def donpclustering(mlnfile,evidfile,E1,E2,E3,outmfile,outefile,outclog):
	basefile = os.path.basename(mlnfile)
	r = random.randint(1,1000)
	generateFeatures(mlnfile,1,1,r)
	tmpspecfile  = str(r)+"-tmp-"+basefile[:basefile.rfind(".")]+".spec"
	ep1 = float(E1)
	ep2 = float(E2)
	ep3 = float(E3)
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
	print(tables)
	print(tableids)        
	print(reducedsizes)

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
	print(pvec)
	print(pindvec)
	print(sizevec)

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
			predtables[tidx][i][idx] = predtables[tidx][i][idx]+1.0
	ifile.close()


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
	#DP-Means
	#while True:
	allclusters = []
	lamdax = 10000
	start_time = time.time()
	maxtime = 1000
	
	while True:
		curr_time = time.time()
		if curr_time - start_time > maxtime:
			break
		Z = []
		allclusters = []
		numclusters = []
		origsizes = []
		for i,t in enumerate(tables):
			idx = prednames.index(t[0])
			sz = len(predtables[idx][tableids[i][0]])
			tmp=[]
			for j in range(0,sz,1):
				tmp.append(0)
			Z.append(tmp)
			numclusters.append(1)
			origsizes.append(1)
		for i,t in enumerate(tables):
			idx = prednames.index(t[0])
			sz = len(predtables[idx][tableids[i][0]])
			origsizes[i] = sz
			numfeatures = len(t)
			mu = []
			tmp = []
			for j in range(0,numfeatures,1):
				tmp.append(0)
			mu.append(tmp)
			for j in range(0,sz,1):
				for k in range(0,numfeatures,1):
					ix = prednames.index(t[k])
					mu[0][k] = mu[0][k] + predtables[ix][tableids[i][k]][j]
			for j in range(0,len(mu[0]),1):
				mu[0][j] = float(mu[0][j])/sz
			#Z= []
			for j in range(0,sz,1):
				Z[i][j]=0
			k = 1
			#print("begin")
			prevsizes = []
			for ns in numclusters:
				prevsizes.append(ns)
			while True:
				converged = True
				for j in range(0,sz,1):
					min_d=-1
					argmin_d = -1
					for c in range(0,len(mu),1):
						#d_ic = calc_dist(,mu[c])
						d_ic = 0
						for k1 in range(0,numfeatures,1):
							ix = prednames.index(t[k1])
							d_ic = d_ic + ((mu[c][k1] - predtables[ix][tableids[i][k1]][j]) ** 2)
						#print(str(mu[c][0])+" "+str(mu[c][1])+" "+str(numfeatures)+" "+str(d_ic))  
						if min_d==-1 or d_ic < min_d:
							min_d = d_ic
							argmin_d = c
					if min_d > lamdax:
						Z[i][j] = k
						k = k + 1
						numclusters[i] = k
						tmp = []
						for k1 in range(0,numfeatures,1):
							ix = prednames.index(t[k1])
							tmp.append(predtables[ix][tableids[i][k1]][j])
						mu.append(tmp)
						converged = False
					else:
						if Z[i][j] != argmin_d:
							Z[i][j] = argmin_d
							converged = False
				if converged:
					break
			'''
			ofile = open("d-"+str(i)+".dat",'w')
			for z in Z:
				ofile.write(str(z)+"\n")
			ofile.close()
			'''
		for i in range(0,len(tableids),1):
			clustervec = []
			print("numlusters="+str(numclusters[i])+" origsize="+str(origsizes[i]))
			for j in range(0,numclusters[i],1):
				tmp = []
				clustervec.append(tmp)
			#print(Z[i])
			for j,z in enumerate(Z[i]):
				clustervec[z].append(j)
			allclusters.append(clustervec)
		
		print("Processing EVidence...")
		totalb = 0
		totalatoms = 0
		totalc = 0
		ZA = 0
		alphas = []
		errs = []
		for i,p in enumerate(pvec):
			ix = prednames.index(p)
			totalval = 0
			if len(pindvec[i])==1:
				totalc = totalc + len(allclusters[pindvec[i][0]])
				for j in range(0,len(allclusters[pindvec[i][0]]),1):
					if len(allclusters[pindvec[i][0]][j])==0:
						alphas.append(0)
						errs.append(0)
						continue
					nevd = len(intersect(sevidtable[ix],allclusters[pindvec[i][0]][j]))
					u1 = len(allclusters[pindvec[i][0]][j])
					e1 = nevd/float(u1)
					totalatoms = totalatoms + u1
					#totalb = totalb + (u1-e1)*e1/float(u1)
					totalval = totalval + e1
					alphas.append(1.0-e1)
					ZA = ZA + 1-e1
					errs.append(nevd)
			elif len(pindvec[i])==2:
				totalc = totalc + len(allclusters[pindvec[i][0]])*len(allclusters[pindvec[i][1]])
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
						if total==0:
							alphas.append(0)
							errs.append(0)
							continue
						#totalb = totalb + (total-count)*count/float(total)
						alphas.append(1.0 - (count/float(total)))
						ZA = ZA + 1-(count/float(total))
						errs.append(count)
						totalatoms = totalatoms + total
						totalval = totalval + count/float(total)
		for l,a in enumerate(alphas):
			#if totalval > 0:
			totalb = totalb + (a/ZA)*errs[l]
		B = totalb
		C = totalc
		print("L = "+str(lamdax)+" B = "+str(B)+" C="+str(C)) 
		#if B < ep1 and C < ep2:
		#	break
		if C > ep2*totalatoms:
			break
		if B < ep1*totalatoms and C < ep2*totalatoms and C > ep3*totalatoms:
			break
		#   sys.exit("No solutions possible")
		#print("L = "+str(lamdax)+" B = "+str(totalb)+" C="+str(totalc))
		#if B < ep1:
		#    lamdax=lamdax*0.75
		#else:
		#    lamdax=lamdax*0.1
		lamdax = lamdax*0.5

	#Log the clusters
	ofile = open(outclog,'w')
	for i in range(0,len(allclusters),1):
		for j in range(0,len(allclusters[i]),1):
			ofile.write(str(allclusters[i][j])+" ")
		ofile.write("\n")
	ofile.close()

	print("Writing EVidence...")        
	ofile=open(outefile,'w')
	for i,p in enumerate(pvec):
		ix = prednames.index(p)
		#ofile1 = open(p+".map",'w')
		if len(pindvec[i])==1:
			for j in range(0,len(allclusters[pindvec[i][0]]),1):
				if len(allclusters[pindvec[i][0]][j])==0:
					continue
				pr = len(intersect(sevidtable[ix],allclusters[pindvec[i][0]][j]))/float(len(allclusters[pindvec[i][0]][j]))
				#pr1 = pr
				ofile.write(str(pr)+" "+p+"("+str(j)+")\n")
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
						continue
					pr = count/float(total)
					ofile.write(str(pr)+" "+p+"("+str(j)+","+str(k)+")\n")
					#ofile1.write(str(total)+"\n")
		#ofile1.close()
	ofile.close()


	'''
	#Alchemy format
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
	#print("Starting Gibbs Sampling...")
	#os.system("../Release/cdlearn cmln.txt cevid.txt cquery.txt 100000 1 MAR res.dat dum.dat")
	#os.system("../Release/cdlearn "+outmfile+" "+outefile+" "+ cevid.txt cquery.txt 100000 1 MAR res.dat dum.dat")
	'''
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
#print(sys.argv[1]+" "+sys.argv[2])
#outmfile = "out.mln"
#outefile = "out.evid"
#outclog = "out.clog"

#donpclustering(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4],sys.argv[5],sys.argv[6],sys.argv[7],sys.argv[8])

