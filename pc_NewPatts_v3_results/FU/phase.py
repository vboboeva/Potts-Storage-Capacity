import numpy
import random
#import ipdb
import scipy
import time
import copy
import os
import re
import sys
from numpy import loadtxt

import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
from pylab import *
import mpl_toolkits.mplot3d.axes3d as axes3d
from scipy.interpolate import interp1d
from scipy.interpolate import InterpolatedUnivariateSpline as inter
from scipy import interpolate

from scipy.interpolate import Rbf
from matplotlib import cm

from scipy.interpolate import griddata
from matplotlib import rc

# the axes attributes need to be set before the call to subplot
plt.rc('text', usetex=True)
plt.rc('font', family='serif', weight='bold', size=18)
plt.rc('xtick.major', size=5, pad=7)
plt.rc('ytick.major', size=5, pad=7)
rcParams['text.latex.preamble'] = [r'\usepackage{sfmath} \boldmath']
plt.rc('xtick', labelsize=18)
plt.rc('ytick', labelsize=18)


N=2000
Cm=200
S=5
a=0.1


filename = 'storage_capacity_dataset_0'
data = loadtxt(filename)

p_values = numpy.arange(20,2020,10)/Cm
p_size = numpy.size(p_values)
print(p_size)
U_values = numpy.array((0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.61, 0.62, 0.63, 0.64, 0.65, 0.66, 0.67, 0.68, 0.69, 0.70, 0.80, 0.90))
U_size = numpy.size(U_values)
print(U_size)

A = numpy.zeros((U_size,p_size))

points = numpy.zeros((p_size*U_size,2))
values = numpy.zeros(p_size*U_size)

column = 9 # put 1 for retrieval, 9 for sparsity
if (column == 1):
	savefig = 'phase_retrieval'
elif (column == 9):
	savefig = 'phase_sparsity'

k=0
for i in range(U_size):
	for j in range(p_size):
		points[k,0] = U_values[i]
		points[k,1] = p_values[j]
		values[k] = data[p_size*i+j,column]
		k=k+1
#print(values[0:p_size])

types1 = ['bicubic']# 'quadric', 'bilinear', 'spline16', 'spline36', 'hanning', 'hamming', 'hermite', 'kaiser', 'catrom', 'gaussian', 'bessel', 'mitchell', 'sinc', 'lanczos']

#print(points)

for i in types1:
	for j in range(10,20,10):
		for k in range(10,20,10):
			
			fig = plt.figure(figsize=(5.5,5))
			grid_svalues, grid_pvalues = numpy.meshgrid(numpy.linspace(min(U_values), max(U_values), j), numpy.linspace(min(p_values), max(p_values), k))
			grid = griddata(points, values, (grid_svalues, grid_pvalues), method='linear')
			plt.imshow(grid, cmap=cm.jet, origin="lower", extent=[min(U_values),max(U_values),min(p_values),max(p_values)], interpolation=i, aspect="auto")
			plt.colorbar(cmap=cm)
			plt.clim(0,1)
			plt.axvline(0.5, color='k')
			plt.xlabel('$U$') 
			plt.ylabel(r'$\alpha$')
			plt.ylim(0.1,6)
			#fig.savefig('j%s_k%s_itype%s.png'%(j,k,i), bbox_inches='tight')
			fig.savefig(savefig, bbox_inches='tight')
			fig.savefig(savefig, format='pdf', bbox_inches='tight')
			plt.close()
	
