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
import scipy.cluster.hierarchy as sch
import scipy.spatial.distance as ssd

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

words = ['bear' ,'cat', 'cow' ,'dog', 'horse', 'arm', 'eye', 'foot', 'hand', 'leg',
 'apartment', 'barn', 'church', 'house', 'igloo', 'arch', 'chimney', 'closet',
 'door', 'window', 'coat', 'dress', 'pants', 'shirt', 'skirt', 'bed', 'chair',
 'desk', 'dresser', 'table', 'ant', 'bee', 'beetle', 'butterfly', 'fly', 'bottle',
 'cup', 'glass', 'knife', 'spoon', 'bell', 'key', 'refrigerator', 'telephone',
 'watch', 'chisel', 'hammer', 'pliers', 'saw', 'screwdriver', 'carrot', 'celery',
 'corn', 'lettuce', 'tomato', 'airplane', 'bicycle', 'car', 'train', 'truck']
 
index = loadtxt('clustering_indices.dat', dtype='int')
words = [ words[i] for i in index]

# Now shuffle features across all nouns
data = loadtxt("nouns_features_matrix.dat")
data = data[index,:]

N_feats = 25
N_words = 60
N_shuffle = 100
N_R=24

# Then replicate
nouns_features=numpy.zeros((60,N_R*N_feats))
for i in range(0,N_feats):
	nouns_features[:,i*N_R:(i+1)*N_R] = numpy.array([data[:,i],]*N_R).transpose()

for j in range(N_shuffle):
	for i in range(N_R*N_feats):
		nouns_features[:,i] = numpy.random.permutation(nouns_features[:,i])

# Re-express nouns as feature vectors
for i in range(N_words):
	N = 0
	for j in range(N_R*N_feats):
		N = N + nouns_features[i,j]**2
	nouns_features[i,:] = nouns_features[i,:]/numpy.sqrt(N)
	
# Correlation matrix
C = numpy.zeros((N_words,N_words))

for i in range(N_words):
	for j in range(N_words):
		d1 = 0
		d2 = 0
		for k in range(N_R*N_feats):
			C[i,j] = C[i,j] + nouns_features[i,k]*nouns_features[j,k]
			d1 = d1 + nouns_features[i,k]*nouns_features[i,k]
			d2 = d2 + nouns_features[j,k]*nouns_features[j,k]
		C[i,j] = C[i,j]/(numpy.sqrt(d1)*numpy.sqrt(d2))
		
numpy.savetxt('C_feats_shuffled_across_all_nouns.dat', C, fmt='%.4f')

fig, ax = plt.subplots(1,1) 
cax = ax.matshow(C, cmap=plt.cm.jet, vmin=0, vmax=1)
#fig.colorbar(cax)
ax.set_xticks(numpy.arange(60))
ax.set_yticks(numpy.arange(60))
ax.set_xticklabels(words, rotation='vertical', fontsize=6)
ax.set_yticklabels(words, fontsize=6)
ax.xaxis.set_ticks_position('top')
ax.yaxis.set_ticks_position('left')

fig.savefig('C_feats_shuffled_across_all_nouns.pdf', bbox_inches='tight')
fig.savefig('C_feats_shuffled_across_all_nouns.png', bbox_inches='tight')
