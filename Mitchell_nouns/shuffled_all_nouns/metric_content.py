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

C = loadtxt("C_feats_shuffled_across_all_nouns.dat")

N_words = numpy.size(C[0,:])
g = 3
print("g=",g)

# Pjoint is the joint
norm = 0
for i in range(N_words):
	for j in range(N_words):
		norm=norm+C[i,j]**g
Pjoint = C**g/norm

numpy.savetxt('Pjoint.dat', Pjoint, fmt='%.4f')


# P is the marginal
P=numpy.zeros((N_words))
for i in range(N_words):
	P[i] = numpy.sum(Pjoint[i,:])

I=0	
for i in range(N_words):
	for j in range(N_words):
		I = I + Pjoint[i,j]*numpy.log(Pjoint[i,j]/(P[i]*P[j]))
I=I/numpy.log(2)

print("I=", I)

#metric content
fcor=0
for i in range(N_words):
	fcor= fcor + Pjoint[i,i]
print("fcor=",fcor)	
Im = ( numpy.log(N_words)+ fcor*numpy.log(fcor) + (1-fcor)*numpy.log(1-fcor) - (1-fcor)*numpy.log(N_words-1) )/numpy.log(2)	
IM = ( numpy.log(N_words) + numpy.log(fcor) )/numpy.log(2)	

print("Im=",Im)
print("IM=",IM)
mcontent = (I-Im)/(IM-Im)

print("mcontent=", mcontent)	


