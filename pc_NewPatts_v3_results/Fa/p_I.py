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

def compute_pc(a_values, p_values, num_datasets):
  
  print(len(a_values))
  print(len(p_values))
  p_cc = numpy.zeros((num_datasets,len(a_values),len(p_values)))
  p_c = numpy.zeros((num_datasets,len(a_values)))
  i=0
  j=0
  data2 = numpy.zeros((num_datasets,len(a_values),len(p_values)))
  for k in range(0,num_datasets):
    data = (loadtxt('apf%.1f_f%.2f/storage_capacity_dataset_%s'%(apf, f, k)))
    for i in range(0,len(a_values)):
      data1 = data[(i)*len(p_values):(i+1)*len(p_values),retrieval]
      #print(numpy.size(data1))
      for j in range(0, len(p_values)):
        data2[k,i,j] = data1[j]
    
  for k in range(0,num_datasets):
    for i in range(0,len(a_values)):
      for j in range(0, len(p_values)):      
        if data2[k,i,j]>= frac_crit:      
          p_cc[k,i,j] = 1
        else:
          p_cc[k,i,j] = 0
      if len((numpy.nonzero(p_cc[k,i,:]))[-1]) != 0:
        p_c[k,i] = p_values[((numpy.nonzero(p_cc[k,i,:]))[-1])[-1]]
      if len((numpy.nonzero(p_cc[k,i,:]))[-1]) == 0:
        p_c[k,i] = 0
  print("p_c=", p_c) 
  return p_c
  
def compute_info(a_values, p_values, num_datasets, i):
	for k in range(0,num_datasets):
		data = (loadtxt('apf%.1f_f%.2f/storage_capacity_dataset_%s'%(apf, f, k)))
		i_c = []
		pi_c = []
		pvals = data[(i)*len(p_values):(i+1)*len(p_values),0]
		info = pvals*data[(i)*len(p_values):(i+1)*len(p_values),7]
		frac_retr = data[(i)*len(p_values):(i+1)*len(p_values),retrieval]
		frac_retr_2ndpatt = data[(i)*len(p_values):(i+1)*len(p_values),retrieval2]
	return pvals, info, frac_retr, frac_retr_2ndpatt

#-----------------------------choose, only this section is to be touched------------------------------------------------------#

# put 1 for 0.7, 2 for 0.8, 3 for 0.9,
retrieval = 1
retrieval2 = 4
frac_crit = 0.5

color = 'green'
p_values = numpy.arange(20,1020,10)
a_values = numpy.linspace(0.1,0.9,9)
print(a_values)
num_datasets = 1

#-------------------------------------------plot-------------------------------------------------------------------------------------------------#

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

fig = plt.figure(figsize=(5.5, 5))
#storage capacity
ax= fig.add_subplot(1,1,1) 
p_c = compute_pc(a_values, p_values, num_datasets)
plt.errorbar(a_values, numpy.mean(p_c/Cm, axis = 0), yerr=numpy.std(p_c/Cm, axis= 0), marker= 'o', label=r"$\displaystyle\ %.1f, %.2f $"%(apf, f), color=color)

ax.set_xlabel('$a$', fontsize=18) 
ax.set_ylabel(r'$\alpha_c$', fontsize=18)
#ax.set_xticks(a_values_ticks)
#ax.set_yticks(p_values_ticks)
ax.tick_params(axis='x', labelsize=14)
ax.tick_params(axis='y', labelsize=14)
#plt.xlim(0, 4)
legend = ax.legend(loc='lower center', fontsize=16)
#plt.grid(True)
plt.tight_layout()
plt.savefig('p_c_storage N=%s Cm=%s S=%s U=%s beta=%s.png'%(N, Cm, S, U, beta ))
plt.savefig('p_c_storage N=%s Cm=%s S=%s U=%s beta=%s.pdf'%(N, Cm, S, U, beta ))


fig = plt.figure(figsize=(5.5, 5))
ax = fig.add_subplot(1,1,1)
#information
for j in range(0,len(a_values)):
	pvals, info, frac_retr, frac_retr_2ndpatt = compute_info(a_values, p_values, num_datasets, j)
	ax.plot(pvals/Cm, info/Cm)
ax.set_xlabel('$p/c_m$', fontsize=18) 
ax.set_ylabel('$pI/c_m$', fontsize=18)
ax.tick_params(axis='x', labelsize=14)
ax.tick_params(axis='y', labelsize=14)
#plt.ylim(0, 3)
legend = ax.legend(loc='best', fontsize=16)
plt.tight_layout()
plt.savefig('info N=%s Cm=%s S=%s U=%s beta=%s.png'%(N, Cm, S, U, beta ))
plt.savefig('info N=%s Cm=%s S=%s U=%s beta=%s.pdf'%(N, Cm, S, U, beta ))


fig = plt.figure(figsize=(5.5, 5))
ax2 = fig.add_subplot(1,1,1)
#ax3 = ax2.twinx()
#fraction retrieved
for j in range(0,len(a_values)):
	pvals, info, frac_retr, frac_retr_2ndpatt = compute_info(a_values, p_values, num_datasets, j)
	ax2.plot(pvals/Cm, frac_retr)
	#ax3.plot(pvals/Cm, frac_retr_2ndpatt)
ax2.set_xlabel('$p$', fontsize=18) 
ax2.set_ylabel(r'fraction retrieved', fontsize=18)
#ax3.set_ylabel(r'fraction retrieved another pattern', fontsize=18)
ax2.tick_params(axis='x', labelsize=14)
ax2.tick_params(axis='y', labelsize=14)
#plt.xlim(0, 4)
legend = ax2.legend(loc='best', fontsize=16)
plt.tight_layout()
plt.savefig('frac retr N=%s Cm=%s S=%s U=%s beta=%s.png'%(N, Cm, S, U, beta ))
plt.savefig('frac retr N=%s Cm=%s S=%s U=%s beta=%s.pdf'%(N, Cm, S, U, beta ))

