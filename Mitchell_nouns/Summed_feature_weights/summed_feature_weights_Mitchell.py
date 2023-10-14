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
from scipy.optimize import curve_fit
from pylab import rcParams

# the axes attributes need to be set before the call to subplot
plt.rc('text', usetex=True)
plt.rc('font', family='serif', weight='bold', size=18)
plt.rc('xtick.major', size=5, pad=7)
plt.rc('ytick.major', size=5, pad=7)
rcParams['text.latex.preamble'] = [r'\usepackage{sfmath} \boldmath']
plt.rc('xtick', labelsize=18)
plt.rc('ytick', labelsize=18)

N_feats = 50
N_words = 60

nouns_features = loadtxt("nouns_features_matrix_both.dat") #_both

print(numpy.shape(nouns_features))

# delete those we don't want
features = numpy.array(('see', 'say', 'taste', 'wear', 'open', 'run', 'near', 'eat', 'hear', 'drive', 'ride', 'touch', 'break', 'enter', 'move', 'listen', 'approach', 'fill', 'clean', 'lift', 'rub', 'smell', 'fear', 'push', 'manipulate', 'get', 'make', 'take', 'use', 'find', 'give', 'put', 'keep', 'live', 'sit', 'stand', 'follow', 'stop', 'grow', 'walk', 'die', 'build', 'close', 'stay', 'fall', 'cut', 'reach', 'kill', 'pass', 'hit'))
to_delete = [] #to_delete = numpy.array(('use', 'manipulate','open', 'near', 'follow', 'push', 'grow', 'fear', 'approach', 'taste', 'smell', 'rub', 'hear', 'stay', 'listen', 'make', 'get', 'find', 'fall'))
print(numpy.shape(features))
for i in range(numpy.size(to_delete)):
	nouns_features = numpy.delete(nouns_features, numpy.where(features == to_delete[i]), 1)
	features = numpy.delete(features, numpy.where(features == to_delete[i]), 0)
print(numpy.shape(features))
print(numpy.shape(nouns_features))

N_feats = numpy.size(features)
# normalize
norm = []
for i in range(N_words):
	x=0
	for j in range(N_feats):
		x=x+nouns_features[i,j]**2
	nouns_features[i,:] = nouns_features[i,:]/numpy.sqrt(x)

# compute summed feature weights
feature_weights = numpy.zeros(N_feats)

for i in range(N_feats):
	for j in range(N_words):
		feature_weights[i] = feature_weights[i] + nouns_features[j,i]

argsorted = numpy.argsort(feature_weights)
argsorted = argsorted[::-1]
weights_sorted = feature_weights[argsorted]	

feats_sorted = []
for i in argsorted:
	feats_sorted = numpy.append(feats_sorted, features[i])

#FIT
def func(x, a, b):
	return a * numpy.exp(-b * x)

#Fit for the parameters a, b, c of the function func
xdata = numpy.arange(N_feats)
ydata = weights_sorted[:]
featsdata = feats_sorted[:] 
popt, pcov = curve_fit(func, xdata, ydata)
print(popt)
		
fig = plt.figure(figsize=(15, 5))
#plt.plot(xdata, func(xdata, *popt), 'r--', lw=3)
plt.plot(xdata, ydata, 'bo' )
plt.xlabel(r'\bf{feature k}') 
plt.ylabel(r'\bf{semantic distinctiveness}') 
plt.xticks(xdata, feats_sorted,rotation=70)
plt.tick_params(axis='x')
plt.tick_params(axis='y')
#plt.text(10,1, '$\sim \exp(-%.3f k)$'%popt[1])
#plt.xlim(0,23)
#plt.grid()
plt.yscale('log', nonposy='clip') #
#plt.xscale('log')
#plt.savefig('corr_distrib_Mitchell.pdf', bbox_inches='tight')
plt.savefig('summed_feature_weights_Mitchell_nonnorm.png', bbox_inches='tight')
plt.savefig('summed_feature_weights_Mitchell_nonnorm.pdf', bbox_inches='tight')
