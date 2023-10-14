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

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

data = loadtxt("nouns_features_matrix.dat")

N_feats = 25
N_words = 60

#correlation matrix
C = numpy.zeros((N_words,N_words))

for i in range(N_words):
	for j in range(N_words):
		d1 = 0
		d2 = 0
		for k in range(N_feats):
			C[i,j] = C[i,j] + data[i,k]*data[j,k]
			d1 = d1 + data[i,k]*data[i,k]
			d2 = d2 + data[j,k]*data[j,k]
		C[i,j] = C[i,j]/(numpy.sqrt(d1)*numpy.sqrt(d2))

numpy.savetxt('C.dat', C, fmt='%.5f')

words = ['bear' ,'cat', 'cow' ,'dog', 'horse', 'arm', 'eye', 'foot', 'hand', 'leg',
 'apartment', 'barn', 'church', 'house', 'igloo', 'arch', 'chimney', 'closet',
 'door', 'window', 'coat', 'dress', 'pants', 'shirt', 'skirt', 'bed', 'chair',
 'desk', 'dresser', 'table', 'ant', 'bee', 'beetle', 'butterfly', 'fly', 'bottle',
 'cup', 'glass', 'knife', 'spoon', 'bell', 'key', 'refrigerator', 'telephone',
 'watch', 'chisel', 'hammer', 'pliers', 'saw', 'screwdriver', 'carrot', 'celery',
 'corn', 'lettuce', 'tomato', 'airplane', 'bicycle', 'car', 'train', 'truck']

fig, ax = plt.subplots(1,1) 
cax = ax.matshow(C, cmap=plt.cm.jet)
fig.colorbar(cax)
ax.set_xticks(numpy.arange(60))
ax.set_yticks(numpy.arange(60))

ax.xaxis.set_ticks_position('top')
ax.yaxis.set_ticks_position('left')

ax.set_xticklabels(words, rotation='vertical', fontsize=6)
ax.set_yticklabels(words, fontsize=6)

fig.savefig('corr_mat_Mitchell.pdf', bbox_inches='tight')
fig.savefig('corr_mat_Mitchell.png', bbox_inches='tight')
