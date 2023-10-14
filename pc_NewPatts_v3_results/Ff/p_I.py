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

f=0.05
Num_fact = 150
apf = 0.4
zeta=0.000001

#---------------------------------------function to compute p_c----------------------------------------------#

def compute_pc(f_values, p_values, num_datasets):
  
  print(len(f_values))
  print(len(p_values))
  p_cc = numpy.zeros((num_datasets,len(f_values),len(p_values)))
  p_c = numpy.zeros((num_datasets,len(f_values)))
  i=0
  j=0
  data2 = numpy.zeros((num_datasets,len(f_values),len(p_values)))
  for k in range(0,num_datasets):
    data = (loadtxt('a%.1f/storage_capacity_dataset_0'%a))
    for i in range(0,len(f_values)):
      data1 = data[(i)*len(p_values):(i+1)*len(p_values),retrieval]
      #print(numpy.size(data1))
      for j in range(0, len(p_values)):
        data2[k,i,j] = data1[j]
    
  for k in range(0,num_datasets):
    for i in range(0,len(f_values)):
      for j in range(0, len(p_values)):      
        if data2[k,i,j]>= frac_crit:      
          p_cc[k,i,j] = 1
        else:
          p_cc[k,i,j] = 0
      if len((numpy.nonzero(p_cc[k,i,:]))[-1]) != 0:
        p_c[k,i] = p_values[((numpy.nonzero(p_cc[k,i,:]))[-1])[-1]]
      if len((numpy.nonzero(p_cc[k,i,:]))[-1]) == 0:
        p_c[k,i] = 0
  #print("p_c=", p_c) 
  return p_c
  
def compute_info(f_values, p_values, num_datasets, i):
	for k in range(0,num_datasets):
		data = (loadtxt('a%.1f/storage_capacity_dataset_0'))
		i_c = []
		pi_c = []
		pvals = data[(i)*len(p_values):(i+1)*len(p_values),0]
		info = pvals*data[(i)*len(p_values):(i+1)*len(p_values),7]
		frac_retr = data[(i)*len(p_values):(i+1)*len(p_values),retrieval]
		frac_retr_2ndpatt = data[(i)*len(p_values):(i+1)*len(p_values),retrieval2]
	return pvals, info, frac_retr, frac_retr_2ndpatt

#-----------------------------choose, only this section is to be touched------------------------------------------------------#

# put 1 for 0.7, 2 for 0.8, 3 for 0.9, 10 for m_avg
retrieval = 1
retrieval2 = 4
frac_crit = 0.5

color = 'green'

f_values = numpy.array((0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10, 0.15, 0.2, 0.4, 0.6, 0.8, 1.0))

num_datasets = 1

#-------------------------------------------plot-------------------------------------------------------------------------------------------------#

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

fig = plt.figure(figsize=(5.5, 5))
#storage capacity
ax= fig.add_subplot(1,1,1) 

a_values = [0.1, 0.3]
p_values = [numpy.arange(20,2020,10), numpy.arange(10,1010,10)]
i=0
for a in a_values:
	print(i)
	p_c = compute_pc(f_values, p_values[:][i], num_datasets)
	plt.errorbar(f_values, numpy.mean(p_c/Cm, axis = 0), yerr=numpy.std(p_c/Cm, axis= 0), marker= 'o', color=color,  label='$a=%.1f$'%a)
	fvals = (loadtxt('a%.1f/corrs_ff_a%.1f.txt'%(a,a)))[:,0]
	pcvals = (loadtxt('a%.1f/corrs_ff_a%.1f.txt'%(a,a)))[:,2]
	plt.plot(fvals, pcvals, '-.', color=color)
	
	i=i+1

plt.ylim(0, 9)
ax.set_xlabel('$f$') 
ax.set_ylabel(r'$\alpha_c$')

plt.savefig('p_c_ff N=%s Cm=%s S=%s U=%s beta=%s.png'%(N, Cm, S, U, beta ), bbox_inches='tight')
plt.savefig('p_c_ff N=%s Cm=%s S=%s U=%s beta=%s.pdf'%(N, Cm, S, U, beta ), bbox_inches='tight')

