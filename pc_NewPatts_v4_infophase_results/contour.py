import numpy
from numpy import loadtxt
import scipy
import time
import copy
import os
import re
import sys
import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib import cm
from matplotlib import rc

from scipy.interpolate import interp1d
from scipy.interpolate import InterpolatedUnivariateSpline as inter
from scipy import interpolate
from scipy.interpolate import Rbf
from scipy.interpolate import griddata
from scipy.ndimage.filters import gaussian_filter

from pylab import rcParams

# the axes attributes need to be set before the call to subplot
plt.rc('text', usetex=True)
plt.rc('font', family='serif', weight='bold', size=18)
plt.rc('xtick.major', size=5, pad=7)
plt.rc('ytick.major', size=5, pad=7)
rcParams['text.latex.preamble'] = [r'\usepackage{sfmath} \boldmath']
plt.rc('xtick', labelsize=18)
plt.rc('ytick', labelsize=18)

N = 2000
Cm = 200
S = 5
U = 0.5
beta = 200
a=0.1

f=0.05
Num_fact = 150
apf = 0.4

def compute_info(p_values, zeta, i):
	data = (loadtxt('data/storage_capacity_zeta%.3f'%(zeta)))
	#print(numpy.shape(data))
	#print(40*"--")
	i_c = []
	pi_c = []
	pvals = data[(i)*len(p_values):(i+1)*len(p_values),0]
	info = data[(i)*len(p_values):(i+1)*len(p_values),7]
	frac_retr = data[(i)*len(p_values):(i+1)*len(p_values),retrieval]
	frac_retr_2ndpatt = data[(i)*len(p_values):(i+1)*len(p_values),retrieval2]
	sparsity = data[(i)*len(p_values):(i+1)*len(p_values),sp]
	#print(frac_retr)
	return pvals, info, frac_retr, frac_retr_2ndpatt, sparsity

#-----------------------------choose, only this section is to be touched------------------------------------------------------#

# put 1 for 0.7, 2 for 0.8, 3 for 0.9, 4 for 0.7corr, 5 for 0.8corr, 6 for 0.9corr, 7 for info, 9 for sparsity
retrieval = 3
retrieval2 = retrieval+3
sp = 9
frac_crit = 0.5

color = 'green'
p_values = numpy.arange(10,2010,10)

f_values = numpy.array(( 0.010, 0.020, 0.030, 0.040, 0.050, 0.060, 0.070, 0.080, 0.090, 0.100, 0.150, 0.20, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.)) 

# DATA OLD 
#zeta_values = numpy.array(( 0.001, 0.002, 0.005, 0.008, 0.009, 0.010, 0.011, 0.012, 0.015, 0.020, 0.050, 0.100, 0.150, 0.200, 0.250, 0.300, 0.350, 0.400, 0.450, 0.500, 0.550, 0.600, 0.650, 0.700, 0.750, 0.800, 0.850, 0.900, 0.950, 1.000, 1.050, 1.100, 1.150, 1.200, 1.250, 1.300, 1.350, 1.400, 1.450, 1.500, 1.550, 1.6, 1.7, 1.8, 1.9, 2.0, 3.0, 4.0, 5.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0))

#zeta_values = numpy.array(( 0.001, 0.002, 0.005, 0.008, 0.009, 0.010, 0.020, 0.050, 0.100, 0.150, 0.200, 0.300, 0.400, 0.500, 0.600, 0.700, 0.800, 0.900, 1.0, 2.0, 3.0, 4.0, 5.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0))

# DATA NEW
zeta_values = numpy.array((0.001, 0.002, 0.003, 0.004, 0.005, 0.006, 0.007, 0.008, 0.009, 0.010, 0.020, 0.030, 0.040, 0.050, 0.060, 0.070, 0.080, 0.090, 0.100, 0.200, 0.300, 0.400, 0.500, 0.600, 0.700, 0.800, 0.900, 1, 2, 3, 4, 5, 7, 20, 30, 40, 50, 60, 70, 80, 90, 100)) 

xsize=4*1.62
ysize=4

#-------------------------------------------------information

meaninfovals = []
for i in range(0, len(zeta_values)):
	for j in range(0, len(f_values)):
		pvals, info, frac_retr, frac_retr_2ndpatt, sparsity = compute_info(p_values, zeta_values[i], j)
		#print(frac_retr)
		x = numpy.where(frac_retr<0.001)
		indices = numpy.where(info[x] == float('inf'))
		noninf = numpy.delete(info[x], indices)
		#print(numpy.size(noninf))
		if (numpy.size(noninf) != 0):
			meaninfovals = numpy.append(meaninfovals, numpy.mean(noninf))
		else:
			meaninfovals = numpy.append(meaninfovals, 0)
coordinates = [(x,y) for x in zeta_values for y in f_values]
#print(coordinates)
z = meaninfovals/Cm
print(numpy.size(z))

fig = plt.figure(figsize=(xsize, ysize))
# define grid.
xi, yi = numpy.meshgrid(numpy.linspace(min(zeta_values),max(zeta_values),2000),numpy.linspace(min(f_values), max(f_values),500))
# grid the data.
zi = griddata(coordinates, z, (xi,yi), method="linear")
#sigma=0.5
#zi = gaussian_filter(zi, sigma)
# contour the gridded data, plotting dots at the randomly spaced data points.
plt.contour(xi,yi,zi,15,linewidths=0.5)
plt.contourf(xi,yi,zi,15,cmap=plt.cm.jet)
plt.colorbar() # draw colorbar
# plot data points.
#plt.scatter(numpy.repeat(zeta_values, len(f_values)), numpy.tile(f_values, len(zeta_values)), 50, z, edgecolor='black', cmap=cm.jet)
x=[0.01, 0.1, 1.0]
y=[0.2, 0.2, 0.2]
plt.scatter(x, y, s=40, alpha=1, marker='o', c='white', zorder=3)
plt.xscale('log')
plt.yscale('log')

plt.axhline(0.05, color='black')
plt.xlim(min(zeta_values),max(zeta_values))
plt.ylim(min(f_values),max(f_values))
plt.xlabel(r'$\zeta$') 
plt.ylabel(r'$f$')
fig.savefig("contour_resI_zeta_vs_f.png", bbox_inches='tight')
fig.savefig("contour_resI_zeta_vs_f.pdf", bbox_inches='tight')
