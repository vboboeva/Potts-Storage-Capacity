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


lines = open( "raw_data.dat", "r" ).readlines()

#print(numpy.size(lines))

N_feats = 26
N_words = 60

Features = ['see', 'say', 'taste', 'wear', 'open', 'run', 'neared', 'eat', 'hear', 'drive', 'ride', 'touch', 'break', 'enter', 'move', 'listen', 'approach', 'fill', 'clean', 'lift', 'rub', 'smell', 'fear', 'push', 'manipulate', '']

#print(numpy.size(Features))

y=[]
for i in range(N_words):
	#print("i=",i)
	for k in range(N_feats):
		print(Features[k])
		for j in range(N_feats):
			feature = (lines[(i*N_feats)+j].split())[0] 
			if (feature == Features[k]):
				print(feature)
				x = re.findall(r"\d+\.\d*",(lines[(i*N_feats)+j]))
				print(x)
				if (x != []):
					y = numpy.append(y, x)
print(y)
z=[]
for i in range(numpy.size(y)):
	z= numpy.append(z, literal_eval(y[i]))
z = numpy.reshape(z, (60,25))

numpy.savetxt('new_data.dat', z, fmt='%.10f')


concepts = []

for i in range(N_words*N_feats):
	find_for = (lines[i].split())[1]
	if (find_for == 'for'):
		concepts = numpy.append(concepts, (lines[i].split())[2])
print(concepts)
numpy.savetxt('concepts_data.dat', concepts, fmt='%s')
