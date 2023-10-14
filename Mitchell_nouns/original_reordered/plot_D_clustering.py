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
import matplotlib.colors as colors
import scipy.cluster.hierarchy as sch
import scipy.spatial.distance as ssd

plt.rc('text', usetex=True)
plt.rc('font', family='serif')


def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    new_cmap = colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(numpy.linspace(minval, maxval, n)))
    return new_cmap

D = loadtxt("D.dat")
C = loadtxt("C.dat")

words = ['bear' ,'cat', 'cow' ,'dog', 'horse', 'arm', 'eye', 'foot', 'hand', 'leg',
 'apartment', 'barn', 'church', 'house', 'igloo', 'arch', 'chimney', 'closet',
 'door', 'window', 'coat', 'dress', 'pants', 'shirt', 'skirt', 'bed', 'chair',
 'desk', 'dresser', 'table', 'ant', 'bee', 'beetle', 'butterfly', 'fly', 'bottle',
 'cup', 'glass', 'knife', 'spoon', 'bell', 'key', 'refrigerator', 'telephone',
 'watch', 'chisel', 'hammer', 'pliers', 'saw', 'screwdriver', 'carrot', 'celery',
 'corn', 'lettuce', 'tomato', 'airplane', 'bicycle', 'car', 'train', 'truck']
 
# shuffle randomly rows of the matrix to see if final outcome changes

##################
n = 0			 #
shuffle = False  #
##################

if shuffle == False:
	n=0
order = numpy.arange(60)

if shuffle == True:
	for j in range(n):
	  numpy.random.shuffle(order)
	  #print(order)
	  D=D[order,:]
	  D=D[:,order]
	  words = [ words[i] for i in order]
	  C=C[order,:]
	  C=C[:,order]
  
# Compute and plot dendrogram.
fig = pylab.figure()
axdendro = fig.add_axes([0.00,0.1,0.3,0.8])
Y = sch.linkage(D[numpy.triu_indices(60,1)], method='centroid')
Z = sch.dendrogram(Y, orientation='left')
#axdendro.set_xticks([])
axdendro.set_yticks([])

# Plot distance matrix.
axmatrix = fig.add_axes([0.4,0.1,0.7,0.8])
index = Z['leaves']

#print(index)
D = D[index,:]
D = D[:,index]
D=D/(numpy.max(D))

cmap = plt.get_cmap('jet')
new_cmap = truncate_colormap(cmap, 0.0, 0.85)


words = [ words[i] for i in index]
im = axmatrix.matshow(D, aspect='auto', origin='lower', cmap=new_cmap)
axmatrix.xaxis.set_ticks_position('top')
axmatrix.yaxis.set_ticks_position('left')
axmatrix.set_xticks(numpy.arange(60))
axmatrix.set_yticks(numpy.arange(60))
axmatrix.set_xticklabels(words, rotation='vertical', fontsize=6)
axmatrix.set_yticklabels(words, fontsize=6)

# Plot colorbar.
axcolor = fig.add_axes([1.15,0.1,0.02,0.8])
pylab.colorbar(im, cax=axcolor)

# Display and save figure.
fig.savefig('dendrogram_D.png', bbox_inches='tight')
fig.savefig('dendrogram_D.pdf', bbox_inches='tight')

numpy.savetxt('clustering_indices.dat', index, fmt='%i')
