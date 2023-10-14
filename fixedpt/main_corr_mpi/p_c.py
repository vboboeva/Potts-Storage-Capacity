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
from matplotlib import colors
from matplotlib import cm
from matplotlib import rc
from pylab import rcParams
from pylab import *

# the axes attributes need to be set before the call to subplot
plt.rc('text', usetex=True)
plt.rc('font', family='serif', weight='bold', size=18)
plt.rc('xtick.major', size=5, pad=7)
plt.rc('ytick.major', size=5, pad=7)
rcParams['text.latex.preamble'] = [r'\usepackage{sfmath} \boldmath']
plt.rc('xtick', labelsize=18)
plt.rc('ytick', labelsize=18)

# get constants from const_potts C++ file

a = loadtxt("corrs_fa.txt")[0:10,0]
pca1 = loadtxt("corrs_fa.txt")[0:10,3]

plt.plot(a,pca1, 'ro-')

plt.xlabel('$f$') 
plt.ylabel(r'$\alpha_c$')
plt.tick_params(axis='x')
plt.tick_params(axis='y')
legend = plt.legend(loc='best')
#plt.grid(True)
plt.ylim(0,9)

plt.savefig('test.png', bbox_inches='tight')
plt.close()
