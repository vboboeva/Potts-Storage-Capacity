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

D = loadtxt("D.dat")

D = D/(numpy.max(D))

words = ['bear' ,'cat', 'cow' ,'dog', 'horse', 'arm', 'eye', 'foot', 'hand', 'leg',
 'apartment', 'barn', 'church', 'house', 'igloo', 'arch', 'chimney', 'closet',
 'door', 'window', 'coat', 'dress', 'pants', 'shirt', 'skirt', 'bed', 'chair',
 'desk', 'dresser', 'table', 'ant', 'bee', 'beetle', 'butterfly', 'fly', 'bottle',
 'cup', 'glass', 'knife', 'spoon', 'bell', 'key', 'refrigerator', 'telephone',
 'watch', 'chisel', 'hammer', 'pliers', 'saw', 'screwdriver', 'carrot', 'celery',
 'corn', 'lettuce', 'tomato', 'airplane', 'bicycle', 'car', 'train', 'truck']

fig, ax = plt.subplots(1,1) 
cax = ax.matshow(D, cmap=plt.cm.jet)
fig.colorbar(cax)
ax.set_xticks(numpy.arange(60))
ax.set_yticks(numpy.arange(60))

ax.xaxis.set_ticks_position('top')
ax.yaxis.set_ticks_position('left')

ax.set_xticklabels(words, rotation='vertical', fontsize=6)
ax.set_yticklabels(words, fontsize=6)

fig.savefig('dist_mat_Mitchell.pdf', bbox_inches='tight')
fig.savefig('dist_mat_Mitchell.png', bbox_inches='tight')
