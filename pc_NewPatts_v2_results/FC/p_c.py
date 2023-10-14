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
lines = open( "const_potts.h", "r" ).readlines()

N = int( (lines[0].split())[2] )
S = int( (lines[1].split())[2] )
a = float( (lines[2].split())[2] )
p_fact = int( (lines[3].split())[2] )
Num_fact = int( (lines[4].split())[2] )
apf = float( (lines[5].split())[2] )
U = float( (lines[6].split())[2] )
beta = int( (lines[7].split())[2] )

print("N=",N)
print("S=",S)
print("a=",a)

#---------------------------------------function to compute p_c----------------------------------------------#

def compute_pc(a_values, p_values, num_datasets, apf, prob):
  
  p_cc = numpy.zeros((num_datasets,len(a_values),len(p_values)))
  p_c = numpy.zeros((num_datasets,len(a_values)))
  i=0
  j=0
  data2 = numpy.zeros((num_datasets,len(a_values),len(p_values)))
  for k in range(0,num_datasets):
    data = (loadtxt('apf%.1f_prob%.2f/storage_capacity_dataset_%s'%(apf,prob, k)))
    for i in range(0,len(a_values)):
      data1 = data[(i)*len(p_values):(i+1)*len(p_values),retrieval]
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
  return p_c

#-----------------------------choose, only this section is to be touched------------------------------------------------------#

# put 1 for 0.7, 2 for 0.8, 3 for 0.9
frac_crit = 0.5
retrieval = 1
C_values = numpy.array((200, 400, 600, 800, 1000, 1200, 1400, 1600, 1900))

apf = [0.4, 1.0, 0.4, 1.0]
prob = [0.05, 0.05, 0.2, 0.2]
colors = ['green', 'blue', 'red' , 'purple']

#apf = [0.4, 0.4, 1.0, 1.0]
#prob = [0.05, 0.2, 0.05, 0.2]
#colors = ['green', 'blue', 'red', 'yellow' ]
p_values = numpy.array((numpy.arange(10,3010,10),numpy.arange(20,1020,10),numpy.arange(20,1020,10),numpy.arange(20,1020,10)))
num_datasets = [1,1,1,1]
C_values_ticks = C_values/N
p_values_ticks = numpy.arange(-20,110,10)
#-------------------------------------------plot-------------------------------------------------------------------------------------------------#

fig = plt.figure(figsize=(5.5, 5))
ax= fig.add_subplot(1,1,1) 
#ax.set_title(r"$\displaystyle\ N=%s, Cm=%s, S=%s, U=%s, \beta=%s$"%(N, Cm, S, U, beta ), fontsize=13)

for i in range(numpy.size(apf)):
  print(apf[i])
  p_c = compute_pc(C_values, p_values[i][:], num_datasets[i], apf[i], prob[i])
  print("p_c=", p_c) 
  if apf[i] == 0.4 and prob[i] == 0.05:
    plt.errorbar(C_values/N, numpy.mean(p_c/C_values, axis = 0), yerr=numpy.std(p_c/C_values, axis= 0), marker= 'o', label=r"$\displaystyle\ a_{p}=%.1f, f=%.2f $"%(apf[i], prob[i]), color=colors[i])
  else:
    plt.errorbar(C_values/N, numpy.mean(p_c/C_values, axis = 0), yerr=numpy.std(p_c/C_values, axis= 0), marker= 'o', label=r"$\displaystyle\ %.1f, %.2f $"%(apf[i], prob[i]), color=colors[i])
plt.axvline(0.1, ls='--', color='black')  
legend = ax.legend(loc='best')
ax.set_xlabel('$c_m/N$') 
ax.set_ylabel(r'$\alpha_c$')
ax.set_xlim(0,1)
ax.set_ylim(0,10)
plt.savefig('p_c_corr_fC N=%s S=%s a=%.2f U=%s beta=%s.png'%(N, S, a, U, beta ), bbox_inches='tight')
plt.savefig('p_c_corr_fC N=%s S=%s a=%.2f U=%s beta=%s.pdf'%(N, S, a, U, beta ), bbox_inches='tight')
plt.close()
