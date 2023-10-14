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

N = 2000
Cm = 200
S = 5
U = 0.5
beta = 200
a=0.1

f=0.05
Num_fact = 150
apf = 0.4

#-----------------------------choose, only this section is to be touched------------------------------------------------------#

# put 1 for 0.7, 2 for 0.8, 3 for 0.9, 4 for 0.7corr, 5 for 0.8corr, 6 for 0.9corr, 7 for info, 9 for sparsity
retrieval = 3
retrieval2 = retrieval+3
sp = 9
frac_crit = 0.5

color = 'green'
p_values = numpy.arange(10,2010,10)
#zeta_values = numpy.array((0.0001, 0.001, 0.01, 0.011, 0.012, 0.013, 0.014, 0.015, 0.016, 0.017, 0.018, 0.019, 0.02, 0.05, 0.1, 0.2, 0.3, 0.4))
zeta_values = numpy.array(( 0.001000, 0.002000, 0.005000, 0.008000, 0.009000, 0.010000, 0.011000, 0.012000, 0.015000, 0.020000, 0.050000, 0.100000, 0.150000, 0.200000, 0.500000))
plot_vals = numpy.array((0.001000, 0.01, 0.02, 0.1))
num_datasets = 1

xsize=4*1.62
ysize=4

#-------------------------------------------------information

fig = plt.figure(figsize=(xsize, ysize))
ax = fig.add_subplot(1,1,1)


data = (loadtxt('apf%.1f_f%.2f/storage_capacity_dataset_%s'%(apf, f, 0)))

fcor = []
I = []
pp = []

k=8
pp = data[k*numpy.size(p_values):(k+1)*numpy.size(p_values),0]
fcor = data[k*numpy.size(p_values):(k+1)*numpy.size(p_values),retrieval]
I = data[k*numpy.size(p_values):(k+1)*numpy.size(p_values),7]

Im = ( numpy.log(pp) + fcor*numpy.log(fcor) + (1-fcor)*numpy.log(1-fcor) - (1-fcor)*numpy.log(pp-1) )/numpy.log(2)	
IM = ( numpy.log(pp) + numpy.log(fcor) )/numpy.log(2)

print(pp)
print(fcor)
print(I)
print(Im)
print(IM)

plt.plot(fcor, I, 'green')
plt.plot(fcor, Im, 'blue')
plt.plot(fcor, IM, 'red')
ax.set_xlabel(r'$\alpha$') 
ax.set_ylabel(r'$I$')
plt.show()



legend = ax.legend(loc='best')
plt.legend()
plt.savefig('mc N=%s Cm=%s S=%s U=%s beta=%s.png'%(N, Cm, S, U, beta ),bbox_inches='tight')
plt.savefig('mc N=%s Cm=%s S=%s U=%s beta=%s.pdf'%(N, Cm, S, U, beta ),bbox_inches='tight')

