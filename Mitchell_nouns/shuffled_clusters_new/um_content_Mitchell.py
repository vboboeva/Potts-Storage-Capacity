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
from matplotlib.colors import LogNorm
import matplotlib as mpl

# the axes attributes need to be set before the call to subplot
plt.rc('text', usetex=True)
plt.rc('font', family='serif', weight='bold', size=18)
plt.rc('xtick.major', size=5, pad=7)
plt.rc('ytick.major', size=5, pad=7)
rcParams['text.latex.preamble'] = [r'\usepackage{sfmath} \boldmath']
plt.rc('xtick', labelsize=18)
plt.rc('ytick', labelsize=18)

#filename = "_feats_shuffled_across_all_nouns"
#filename = "_feats_shuffled_within_clusters_offblock_mean"


#filename = "_feats_shuffled_within_clusters"
#filename = "_feats_shuffled_within_clusters_offblock_mean_inblock_mean"
filename = '_feats_shuffled_within_clusters_intercluster_mean'
#filename = "_original"

D = loadtxt("D%s.dat"%filename)
N_words = 60
power = 1
#triplets

delta1 = []
delta2 = []

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

numpy.set_printoptions(threshold=numpy.nan)	
print("delta1=",numpy.size(delta1))
print("delta2=",numpy.size(delta2))

# ------------------------------------- compute UM content
lamda = []

for i in range(numpy.size(delta1)):
	if( (numpy.log(delta1[i]) + numpy.log(delta2[i])) != 0 ):
		x = (numpy.log(delta1[i]) - numpy.log(delta2[i]))/(numpy.log(delta1[i]) + numpy.log(delta2[i]))
		lamda = numpy.append(lamda,x)

fig = plt.figure(figsize=(5, 5))
plt.hist(lamda, normed=1, bins='auto')
plt.axvline(numpy.mean(lamda))
plt.suptitle(r'$\lambda=%.2f$'%(numpy.mean(lamda)))
plt.savefig('lamda_distrib%s_dpow%.1f.png'%(filename,power), bbox_inches='tight')
print("umcontent=", numpy.mean(lamda))

# ------------------------------------- PLOT

fig = plt.figure(figsize=(3.5, 5))
#plt.plot(numpy.linspace(0.0,1.1,11),numpy.linspace(0.0,1.1,11), alpha=1, color='black', lw=3, zorder=2)
#plt.plot(numpy.linspace(0.0,1.1,11),1-numpy.linspace(0.0,1.1,11), alpha=1, color='black', lw=2, zorder=2)
#plt.axvline(1.0, alpha=1, color='black', zorder=3)
#plt.scatter(delta2, delta1, marker='o',s=0.1, c='black', zorder=3)
#plt.hexbin(delta2, delta1, extent=[0.5, 1.1, 0.0, 1.1], gridsize=25, bins='log', cmap='jet')
#####
x = numpy.linspace(0.0,1.1,11)
y = x**(1+numpy.mean(lamda)/(1-numpy.mean(lamda)))
plt.plot(x,y,'r--')
#plt.axvline(1.0, color='red', ls='--', zorder=3)

plt.plot(numpy.linspace(0.0,1.1,11),numpy.linspace(0.0,1.1,11), alpha=0.5, color='black', lw=1, zorder=2)
plt.plot(numpy.linspace(0.0,1.1,11),1-numpy.linspace(0.0,1.1,11), alpha=0.5, color='black', lw=1, zorder=2)
plt.hist2d(delta2, delta1, range=numpy.array([(0.5, 1.1), (0, 1.1)]), bins=30, norm=LogNorm(),cmap='Greys', zorder=1)
plt.xlim(0.5, 1.1)
plt.ylim(0, 1)
plt.clim(1,1100)
plt.colorbar()
#####
plt.xlabel('$\delta_2 (d_{med}/d_{max})$') 
plt.ylabel('$\delta_1 (d_{min}/d_{max})$')

plt.savefig('um_content%s_dpow%.1f.pdf'%(filename,power), bbox_inches='tight')
plt.savefig('um_content%s_dpow%.1f.png'%(filename,power), bbox_inches='tight')

