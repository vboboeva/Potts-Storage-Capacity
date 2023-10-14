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

alphavals=[4,5,6,6.30,6.40,6.50]

for alpha in alphavals:
	time = loadtxt("test_a%.2f_alpha%.2f.txt"%(a,alpha))[:,0]
	m = loadtxt("test_a%.2f_alpha%.2f.txt"%(a,alpha))[:,1]
	plt.plot(time,m, label="%.2f"%alpha)

plt.axhline(1)

ax.tick_params(axis='x')
ax.tick_params(axis='y')
ax.set_xlabel(r"$a$")
ax.set_ylabel(r"$\alpha_c$")  
ax.set_ylim(0,1.1)
plt.legend(loc='best', frameon=False)
plt.savefig('alpha_Fa.png', bbox_inches='tight')
plt.show()
