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

nouns_features = loadtxt("nouns_features_matrix.dat")

N_feats = 25
N_words = 60

#correlation matrix
C = numpy.zeros((N_words,N_words))

for i in range(N_words):
	for j in range(N_words):
		d1 = 0
		d2 = 0
		for k in range(N_feats):
			C[i,j] = C[i,j] + nouns_features[i,k]*nouns_features[j,k]
			d1 = d1 + nouns_features[i,k]*nouns_features[i,k]
			d2 = d2 + nouns_features[j,k]*nouns_features[j,k]
		C[i,j] = C[i,j]/(numpy.sqrt(d1)*numpy.sqrt(d2))
			
fig = plt.figure(figsize=(5, 5))
plt.hist(C[numpy.triu_indices(60, 1)], bins='auto', normed=1)
plt.xlabel('$correlation$') 
plt.ylabel('$pdf$') 
plt.tick_params(axis='x')
plt.tick_params(axis='y')
plt.savefig('corr_distrib_Mitchell.pdf', bbox_inches='tight')
plt.savefig('corr_distrib_Mitchell.png', bbox_inches='tight')
