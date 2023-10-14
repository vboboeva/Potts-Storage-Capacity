# -*- coding: utf-8 -*-
import numpy
from numpy import loadtxt
import scipy
import pylab
import time
import copy
import os
import re
import sys
import matplotlib.pyplot as plt
from matplotlib import rc
from scipy import stats
from itertools import combinations
from ast import literal_eval
from matplotlib.colors import ListedColormap, NoNorm
from pylab import rcParams

# the axes attributes need to be set before the call to subplot
plt.rc('text', usetex=True)
plt.rc('font', family='serif', weight='bold', size=18)
plt.rc('xtick.major', size=5, pad=7)
plt.rc('ytick.major', size=5, pad=7)
rcParams['text.latex.preamble'] = [r'\usepackage{sfmath} \boldmath']
plt.rc('xtick', labelsize=18)
plt.rc('ytick', labelsize=18)

#filename = "_feats_shuffled_across_all_nouns.dat"
#filename = "_feats_shuffled_within_clusters.dat"
#filename = "_feats_shuffled_within_clusters_offblock_mean.dat"
filename = "_feats_shuffled_within_clusters_offblock_mean_inblock_mean.dat"
#filename = ".dat"

D = loadtxt("D%s"%filename)
N_words = 60
power = 1
#triplets

delta1 = []
delta2 = []

i_else = numpy.arange(5,11+1)

D=numpy.delete(D, numpy.s_[5:11+1], 0)
D=numpy.delete(D, numpy.s_[5:11+1], 1)

N_words = numpy.shape(D)[0]
for item in list(combinations(range(N_words),3)):
	#print(item)
	d1 = D[item[0],item[1]]#numpy.fabs()
	d2 = D[item[1],item[2]]#numpy.fabs()
	d3 = D[item[2],item[0]]#numpy.fabs()
	
	#print(d1)
	#print(d2)
	#print(d3)
	dmin = (numpy.sort([d1,d2,d3]))[0]
	dint = (numpy.sort([d1,d2,d3]))[1]
	dmax = (numpy.sort([d1,d2,d3]))[2]
	#print("dmin=",dmin)
	#print("dint=",dint)
	#print("dmax=",dmax)
	#print(10*"--")
	
	if (dmin > 0 and dint > 0 and dmax > 0):
		delta1 = numpy.append(delta1, dmin/dmax)
		delta2 = numpy.append(delta2, dint/dmax)
	
#print("delta1=",delta1)
#print("delta2=",delta2)

# ------------------------------------- compute UM content
lamda = []

for i in range(numpy.size(delta1)):
	if( (numpy.log(delta1[i]) + numpy.log(delta2[i])) != 0 ):
		x = (numpy.log(delta1[i]) - numpy.log(delta2[i]))/(numpy.log(delta1[i]) + numpy.log(delta2[i]))
		lamda = numpy.append(lamda,x)
print(lamda)
print("umcontent=", numpy.mean(lamda))

# ------------------------------------- PLOT

fig = plt.figure(figsize=(3.5, 5))
plt.plot(numpy.linspace(0.0,1.1,11),numpy.linspace(0.0,1.1,11), alpha=1, color='white', lw=3, zorder=2)
plt.plot(numpy.linspace(0.0,1.1,11),1-numpy.linspace(0.0,1.1,11), alpha=1, color='white', lw=2, zorder=2)
#plt.axvline(1.0, alpha=1, color='black', zorder=3)
plt.scatter(delta2, delta1, marker='o',s=0.1, c='black', zorder=3)
#####
#plt.hexbin(numpy.append(numpy.linspace(0.0,0.8),delta2), numpy.append(numpy.linspace(0.0,0.8), delta1), gridsize=50, bins='log', cmap='jet', vmin=0, vmax=2)
#plt.colorbar()
#####
plt.xlabel('$\delta_2 (d_{med}/d_{max})$') 
plt.ylabel('$\delta_1 (d_{min}/d_{max})$')
plt.xlim(0.5,1.1)
plt.ylim(0.,1.1)

print(list(combinations(range(8),2)))

#plt.savefig('um_content%s_dpow%.1f.pdf'%(filename,power), bbox_inches='tight')
#plt.savefig('um_content%s_dpow%.1f.png'%(filename,power), bbox_inches='tight')

