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

D = loadtxt("D.dat")
C = loadtxt("C.dat")

words = ['bear' ,'cat', 'cow' ,'dog', 'horse', 'arm', 'eye', 'foot', 'hand', 'leg',
 'apartment', 'barn', 'church', 'house', 'igloo', 'arch', 'chimney', 'closet',
 'door', 'window', 'coat', 'dress', 'pants', 'shirt', 'skirt', 'bed', 'chair',
 'desk', 'dresser', 'table', 'ant', 'bee', 'beetle', 'butterfly', 'fly', 'bottle',
 'cup', 'glass', 'knife', 'spoon', 'bell', 'key', 'refrigerator', 'telephone',
 'watch', 'chisel', 'hammer', 'pliers', 'saw', 'screwdriver', 'carrot', 'celery',
 'corn', 'lettuce', 'tomato', 'airplane', 'bicycle', 'car', 'train', 'truck']

index = loadtxt('clustering_indices.dat', dtype='int')
words = [ words[i] for i in index]
# Now plot correlation matrix according to ordering given by the clustering
fig, ax = plt.subplots(1,1) 
C = C[index,:]
C = C[:,index]
cax = ax.matshow(C, cmap=plt.cm.jet, vmin=0, vmax=1)
#fig.colorbar(cax)
ax.set_xticks(numpy.arange(60))
ax.set_yticks(numpy.arange(60))
ax.xaxis.set_ticks_position('top')
ax.yaxis.set_ticks_position('left')
ax.set_xticklabels(words, rotation='vertical', fontsize=6)
ax.set_yticklabels(words, fontsize=6)

# Display and save figure.
fig.savefig('C_original.png', bbox_inches='tight')
fig.savefig('C_original.pdf', bbox_inches='tight')

# save reordered C
numpy.savetxt('C_original_reordered.dat', C, fmt='%.5f')

