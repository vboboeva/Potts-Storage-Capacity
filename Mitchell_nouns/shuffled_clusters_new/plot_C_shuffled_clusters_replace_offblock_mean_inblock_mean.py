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
from itertools import permutations
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

C = loadtxt('C_feats_shuffled_within_clusters.dat')

# Define cluster indices
i_truck_arm = [53,54,55,56,57,58,59]
i_tomato_eye = [39,40,41,42,43,44,45,46,47,48,49,50,51,52]
i_lettuce_cow = [32,33,34,35,36,37,38]
i_refr_screw = [19,20,21,22,23,24,25,26,27,28,29,30,31]
i_car_beetle = [16,17,18]
i_pants_shirt = [12,13,14,15] 
i_glass = [11]
i_watch = [10]
i_key = [9]
i_apartment = [8]
i_ant = [7]
i_chimney = [6]
i_hammer = [5]
i_bike_train = [1,2,3,4]
i_carrot = [0]

# compute mean of off-block entries

num_on_block=0

ivals = []
jvals = []

mean_truck_arm, mean_tomato_eye, mean_lettuce_cow, mean_refr_screw, mean_car_beetle, mean_bike_train, mean_pants_shirt= [],[],[],[],[],[],[]

for i in range(N_words):
	for j in range(N_words):
		if (i in i_truck_arm and j in i_truck_arm):
			num_on_block=num_on_block+1
			ivals = numpy.append(ivals,i) 
			jvals = numpy.append(jvals,j)
			if (i!=j):
				mean_truck_arm = numpy.append(mean_truck_arm, C[i,j])
			
		elif (i  in i_tomato_eye and j  in i_tomato_eye):
			num_on_block=num_on_block+1
			ivals = numpy.append(ivals,i) 
			jvals = numpy.append(jvals,j) 
			if (i!=j):			
				mean_tomato_eye = numpy.append(mean_tomato_eye, C[i,j])
			
		elif (i  in i_lettuce_cow and j  in i_lettuce_cow):
			num_on_block=num_on_block+1
			ivals = numpy.append(ivals,i) 
			jvals = numpy.append(jvals,j)
			if (i!=j):	
				mean_lettuce_cow = numpy.append(mean_lettuce_cow, C[i,j])
		
		elif (i  in i_refr_screw and j  in i_refr_screw):
			num_on_block=num_on_block+1
			ivals = numpy.append(ivals,i) 
			jvals = numpy.append(jvals,j) 		
			if (i!=j):
				mean_refr_screw = numpy.append(mean_refr_screw, C[i,j])
			
		elif (i  in i_car_beetle and j  in i_car_beetle):
			num_on_block=num_on_block+1
			ivals = numpy.append(ivals,i) 
			jvals = numpy.append(jvals,j) 		
			if (i!=j):
				mean_car_beetle = numpy.append(mean_car_beetle, C[i,j])
			
		elif (i  in i_bike_train and j  in	i_bike_train): 
			num_on_block=num_on_block+1
			ivals = numpy.append(ivals,i) 
			jvals = numpy.append(jvals,j)
			if (i!=j): 		
				mean_bike_train = numpy.append(mean_bike_train, C[i,j])
			
		elif (i in i_pants_shirt and j in i_pants_shirt): 
			num_on_block=num_on_block+1
			ivals = numpy.append(ivals,i) 
			jvals = numpy.append(jvals,j) 			
			if (i!=j):
				mean_pants_shirt = numpy.append(mean_pants_shirt, C[i,j])
			
# set off-block to mean
print(num_on_block)
on_block_couples = list(zip(ivals,jvals))
all_couples = permutations(numpy.arange(N_words), 2)

a=[]
for combo in all_couples:
	#print(combo)
	if combo not in on_block_couples:
		a = numpy.append(a, C[combo[0], combo[1]])
mymean = numpy.mean(a)
print(mymean)

on_block_couples = list(zip(ivals,jvals))
all_couples = permutations(numpy.arange(N_words), 2)

for combo1 in all_couples:
	#print(combo1)
	if combo1 not in on_block_couples:
		C[combo1[0], combo1[1]] = mymean
		
# set in-block to mean

for i in range(N_words):
	for j in range(N_words):
		if (i in i_truck_arm and j in i_truck_arm and i!=j):
			C[i,j] = numpy.mean(mean_truck_arm)
			
		elif (i  in i_tomato_eye and j  in i_tomato_eye and i!=j):
			C[i,j] = numpy.mean(mean_tomato_eye)
			
		elif (i  in i_lettuce_cow and j  in i_lettuce_cow and i!=j):
			C[i,j] = numpy.mean(mean_lettuce_cow)
			
		elif (i  in i_refr_screw and j  in i_refr_screw and i!=j):
			C[i,j] = numpy.mean(mean_refr_screw)
			
		elif (i  in i_car_beetle and j  in i_car_beetle and i!=j):
			C[i,j] = numpy.mean(mean_car_beetle)
			
		elif (i  in i_bike_train and j  in	i_bike_train and i!=j): 
			C[i,j] = numpy.mean(mean_bike_train)
			
		elif (i in i_pants_shirt and j in i_pants_shirt and i!=j): 
			C[i,j] = numpy.mean(mean_pants_shirt)
		
	
numpy.savetxt('C_feats_shuffled_within_clusters_offblock_mean_inblock_mean.dat', C, fmt='%.4f')
	
#----------------------plot

fig, ax = plt.subplots(1,1) 
cax = ax.matshow(C, cmap=plt.cm.jet, vmin=0, vmax=1)
fig.colorbar(cax)
ax.set_xticks(numpy.arange(N_words))
ax.set_yticks(numpy.arange(N_words))
ax.set_xticklabels(words, rotation='vertical', fontsize=6)
ax.set_yticklabels(words, fontsize=6)
ax.xaxis.set_ticks_position('top')
ax.yaxis.set_ticks_position('left')

fig.savefig('C_clustered_shuffled_within_%s_offblock_mean_inblock_mean.pdf'%N_shuffle, bbox_inches='tight')
fig.savefig('C_clustered_shuffled_within_%s_offblock_mean_inblock_mean.png'%N_shuffle, bbox_inches='tight')
