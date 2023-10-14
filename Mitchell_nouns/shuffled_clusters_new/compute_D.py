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

#filename = '_feats_shuffled_across_all_nouns'
#filename = '_feats_shuffled_within_clusters'
#filename = '_feats_shuffled_within_clusters_offblock_mean'
filename = '_feats_shuffled_within_clusters_offblock_mean_inblock_mean'
#filename = '_feats_shuffled_within_clusters_intercluster_mean'
#filename = ''

C = loadtxt("C%s.dat"%filename)

power = 0.3

N_feats = 25
N_words = 60

PP= numpy.zeros((N_words,N_words))
for i in range(N_words):
	norm = numpy.sum(C[i,:])
	#print(norm)
	for j in range(N_words):
		PP[i,j] = C[i,j]/norm
		
#print(PP)
D= numpy.zeros((N_words,N_words))
for i in range(N_words):
	for j in range(N_words):
		if PP[i,j]!=0 and PP[j,i]!=0 and PP[i,i]!=0 and PP[j,j]!=0:
			D[i,j] = - numpy.log(PP[i,j]*PP[j,i]/(PP[i,i]*PP[j,j]))
			D[i,j] = pow(D[i,j], power)	

print(D)
numpy.savetxt("D%s.dat"%filename, D, fmt='%.4f')
