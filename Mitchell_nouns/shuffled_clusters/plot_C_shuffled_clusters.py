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

# Shuffle the feature to noun weights within each of the below-defined clusters:

#a) truck to arm words [53:59]
#b) tomato to eye [39:52]
#c) lettuce to cow [32:38]
#d) refrigerator to screwdriver [19:31]
#e) two mini-blocks: car to beetle [16:18] plus bicycle to train [1:4]
#f) pants to shirt (the red minicluster) [12:15]
#g) everything else [5:11]
#h) carrot, corn, celery [0] [36:37] 

#print(words.index("celery"))

# Reorder noun-feature matrix according to clustering
N_feats = 25
N_words = 60
N_shuffle = 10
N_R = 100

data = loadtxt("nouns_features_matrix.dat")
data = data[index,:]

# Then replicate

nouns_features=numpy.zeros((N_words,N_R*N_feats))
for i in range(0,N_feats):
	nouns_features[:,i*N_R:(i+1)*N_R] = numpy.array([data[:,i],]*N_R).transpose()

# Define cluster indices
i_truck_arm = numpy.arange(53,59+1)
i_tomato_eye = numpy.arange(39,52+1)
i_lettuce_cow = numpy.arange(32,38+1)
i_refr_screw = numpy.arange(19,31+1)
i_car_beetle = numpy.arange(16,18+1)
i_bike_train = numpy.arange(1,4+1)
i_else = numpy.arange(5,11+1)
i_veggies = numpy.arange(36,37+1)


# Make an average
C = numpy.zeros((N_words,N_words))
# Shuffle
for j in range(N_shuffle):
	
	for i in range(N_R*N_feats):
		nouns_features[i_truck_arm,i] = numpy.random.permutation(nouns_features[i_truck_arm,i])
		nouns_features[i_tomato_eye,i] = numpy.random.permutation(nouns_features[i_tomato_eye,i])
		nouns_features[i_lettuce_cow,i] = numpy.random.permutation(nouns_features[i_lettuce_cow,i])
		nouns_features[i_refr_screw,i] = numpy.random.permutation(nouns_features[i_refr_screw,i])
		nouns_features[i_car_beetle,i] = numpy.random.permutation(nouns_features[i_car_beetle,i])
		nouns_features[i_bike_train,i] = numpy.random.permutation(nouns_features[i_bike_train,i])
		nouns_features[i_else,i] = numpy.random.permutation(nouns_features[i_else,i])
		nouns_features[i_veggies,i] = numpy.random.permutation(nouns_features[i_veggies,i])
		
	# Re-express nouns as feature vectors
	for i in range(N_words):
		N = 0
		for j in range(N_R*N_feats):
			N = N + nouns_features[i,j]**2
		nouns_features[i,:] = nouns_features[i,:]/numpy.sqrt(N)

	C1 = numpy.zeros((N_words,N_words))
	# Recompute correlation matrix
	for i in range(N_words):
		for j in range(N_words):
			d1 = 0
			d2 = 0
			for k in range(N_R*N_feats):
				C1[i,j] = C1[i,j] + nouns_features[i,k]*nouns_features[j,k]
				d1 = d1 + nouns_features[i,k]*nouns_features[i,k]
				d2 = d2 + nouns_features[j,k]*nouns_features[j,k]
			C1[i,j] = C1[i,j]/(numpy.sqrt(d1)*numpy.sqrt(d2))
			C[i,j] = C[i,j]+C1[i,j]
print(C)
C[:,:] = C[:,:]/N_shuffle	
print(C)

#----------------------this is to remove the null cluster for the ultrametric cluster
# C_remove_lowcorr_cluster = C
# for i in i_else:
	# print(i)
	# C_remove_lowcorr_cluster = numpy.delete(C_remove_lowcorr_cluster,i,0)
	# C_remove_lowcorr_cluster = numpy.delete(C_remove_lowcorr_cluster,i,1)

# C_remove_lowcorr_cluster = numpy.delete(C_remove_lowcorr_cluster,0,0)
# C_remove_lowcorr_cluster = numpy.delete(C_remove_lowcorr_cluster,0,1)

# print('this=',numpy.shape(C_remove_lowcorr_cluster))

numpy.savetxt('C_feats_shuffled_within_clusters.dat', C, fmt='%.4f')

#----------------------this is to remove the null cluster

fig, ax = plt.subplots(1,1) 
cax = ax.matshow(C, cmap=plt.cm.jet, vmin=0, vmax=1)
#fig.colorbar(cax)
ax.set_xticks(numpy.arange(N_words))
ax.set_yticks(numpy.arange(N_words))
ax.set_xticklabels(words, rotation='vertical', fontsize=6)
ax.set_yticklabels(words, fontsize=6)
ax.xaxis.set_ticks_position('top')
ax.yaxis.set_ticks_position('left')

fig.savefig('C_clustered_shuffled_within_%s.pdf'%N_shuffle, bbox_inches='tight')
fig.savefig('C_clustered_shuffled_within_%s.png'%N_shuffle, bbox_inches='tight')
