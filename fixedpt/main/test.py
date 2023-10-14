import numpy
import matplotlib.pyplot as plt
import scipy
import copy
import os
import re
import sys
from matplotlib import rc
from numpy import loadtxt
from scipy.stats import binom
from pylab import *
from matplotlib.font_manager import FontProperties

# the axes attributes need to be set before the call to subplot
plt.rc('text', usetex=True)
plt.rc('font', family='serif', weight='bold', size=18)
plt.rc('xtick.major', size=5, pad=7)
plt.rc('ytick.major', size=5, pad=7)
rcParams['text.latex.preamble'] = [r'\usepackage{sfmath} \boldmath']
plt.rc('xtick', labelsize=18)
plt.rc('ytick', labelsize=18)

fig = plt.figure(figsize=(5, 5))
ax= fig.add_subplot(1,1,1) 


a=0.2

alpha=9.1

m=loadtxt("alpha%.2f.txt"%(alpha))[:,1]
print(m)
plt.hist(m, normed=1, bins=20)

ax.tick_params(axis='x')
ax.tick_params(axis='y')
ax.set_xlabel(r"$m$")
#ax.set_ylim(0,1.1)
#plt.legend(loc='best', frameon=False)
plt.savefig('m_distrib_alpha%.2f.png'%alpha, bbox_inches='tight')
