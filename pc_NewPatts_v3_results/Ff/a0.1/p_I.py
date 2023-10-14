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
    data = (loadtxt('storage_capacity_dataset_0'))
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
		data = (loadtxt('storage_capacity_dataset_0'))
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
p_values = numpy.arange(20,2020,10)

f_values = numpy.array((0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10, 0.15, 0.2, 0.4, 0.6, 0.8, 1.0))

print(f_values)
num_datasets = 1

#-------------------------------------------plot-------------------------------------------------------------------------------------------------#

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

fig = plt.figure(figsize=(5.5, 5))
#storage capacity
ax= fig.add_subplot(1,1,1) 
p_c = compute_pc(f_values, p_values, num_datasets)
print(p_c/Cm)
plt.errorbar(f_values, numpy.mean(p_c/Cm, axis = 0), yerr=numpy.std(p_c/Cm, axis= 0), marker= 'o', color=color)
plt.xlim(0, 1)
plt.ylim(0, 8)
ax.set_xlabel('$f$') 
ax.set_ylabel(r'$\alpha_c$')


fvals = (loadtxt('corrs_ff_a0.1.txt'))[:,0]
pcvals = (loadtxt('corrs_ff_a0.1.txt'))[:,2]

plt.plot(fvals, pcvals, '-.', color='green')

plt.savefig('p_c_ff N=%s Cm=%s S=%s U=%s beta=%s.png'%(N, Cm, S, U, beta ), bbox_inches='tight')
plt.savefig('p_c_ff N=%s Cm=%s S=%s U=%s beta=%s.pdf'%(N, Cm, S, U, beta ), bbox_inches='tight')


fig = plt.figure(figsize=(5.5, 5))
ax = fig.add_subplot(1,1,1)
#information
for j in range(0,len(f_values)):
	pvals, info, frac_retr, frac_retr_2ndpatt = compute_info(f_values, p_values, num_datasets, j)
	ax.plot(pvals/Cm, info/Cm)
ax.set_xlabel('$p/c_m$') 
ax.set_ylabel('$pI/c_m$')
#plt.ylim(0, 3)
legend = ax.legend(loc='best')
plt.tight_layout()
#plt.savefig('info N=%s Cm=%s S=%s U=%s beta=%s.png'%(N, Cm, S, U, beta ), bbox_inches='tight')
#plt.savefig('info N=%s Cm=%s S=%s U=%s beta=%s.pdf'%(N, Cm, S, U, beta ), bbox_inches='tight')


fig = plt.figure(figsize=(5.5, 5))
ax2 = fig.add_subplot(1,1,1)
#ax3 = ax2.twinx()
#fraction retrieved
for j in range(0,len(f_values)):
	pvals, info, frac_retr, frac_retr_2ndpatt = compute_info(f_values, p_values, num_datasets, j)
	ax2.plot(pvals/Cm, frac_retr)
	#ax3.plot(pvals/Cm, frac_retr_2ndpatt)
ax2.set_xlabel('$p$') 
ax2.set_ylabel(r'fraction retrieved')
legend = ax2.legend(loc='best')
plt.tight_layout()
#plt.savefig('frac retr N=%s Cm=%s S=%s U=%s beta=%s.png'%(N, Cm, S, U, beta ), bbox_inches='tight')
#plt.savefig('frac retr N=%s Cm=%s S=%s U=%s beta=%s.pdf'%(N, Cm, S, U, beta ), bbox_inches='tight')

