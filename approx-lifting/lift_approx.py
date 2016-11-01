import os
import sys
#import matplotlib
import matplotlib.pyplot as plt
import npclustering as npc
import pclustering as pc
import numpy as np
from collections import defaultdict
import math
import matplotlib.gridspec as gridspec

if len(sys.argv)<5:
	sys.exit("Usage: lift_approx -np/-p mln evidence query [cratio required for -p] [clusterUbound required for -np] [clusterLbound required for -np]")
if sys.argv[1]=="-p":
	pc.doclustering(sys.argv[2],sys.argv[3],sys.argv[5],"cmln.txt","cevid.txt","clog.dat")
else:
	npc.donpclustering(sys.argv[2],sys.argv[3],0.1,sys.argv[5],sys.argv[6],"cmln.txt","cevid.txt","clog.dat")

ifile = open("cevid.txt")
#Y = []
pred_dict = defaultdict(list)
for l in ifile:
	parts = l.strip().split()
	pred_dict[parts[1][:parts[1].find("(")]].append(float(parts[0]))
	#Y.append(float(parts[0]))
ifile.close()
v1 = math.ceil(len(pred_dict.keys())/2)
v2 = math.ceil(len(pred_dict.keys())/v1)
gs = gridspec.GridSpec(v1,v2)
c1 = 0
c2 = 0

for i,k in enumerate(pred_dict.keys()):
	#plt.subplot(c1,c2%v,(i%v)+1)
	plt.subplot2grid((v1,v2),(c1,c2))
	c2 = c2 + 1
	if c2==v2:
		c1 = (c1 + 1)%v1
		c2 = 0
	plt.scatter(np.arange(0,len(pred_dict[k]),1),pred_dict[k],c=pred_dict[k])
	plt.xlabel(k)
	plt.ylabel("Evidence")
	plt.ylim((0,1))
	plt.xlim((0,len(pred_dict[k])))
	plt.colorbar()
#plt.show()
#plt.title(plotname)
#plt.show()

#plt.xlabel("Meta-Atoms")
#plt.ylabel("Evidence Probability")
#plt.scatter(np.arange(0,len(Y),1),Y,c=Y)
#plt.title(plotname)
#plt.colorbar()
#plt.ylim((0,1))
#plt.xlim((0,len(Y)))
plt.tight_layout()
plt.show()
plt.close()