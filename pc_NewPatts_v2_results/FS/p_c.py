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
Cm = int( (lines[1].split())[2] )
a = float( ( lines[2].split())[2] )
U = float( (lines[6].split())[2] )
beta = int( (lines[7].split())[2] )

print("N=",N)
print("Cm=",Cm)
print("a=",a)
print("U=",U)
print("beta=",beta)

#---------------------------------------function to compute p_c----------------------------------------------#

def compute_pc(S_values, p_values, num_datasets, apf, prob):
  
  p_cc = numpy.zeros((num_datasets,len(S_values),len(p_values)))
  p_c = numpy.zeros((num_datasets,len(S_values)))
  i=0
  j=0
  data2 = numpy.zeros((num_datasets,len(S_values),len(p_values)))
  for k in range(0,num_datasets):
    data = (loadtxt('apf%.2f_prob%.2f/storage_capacity_dataset_%s'%(apf, prob, k)))
    for i in range(0,len(S_values)):
      data1 = data[(i)*len(p_values):(i+1)*len(p_values),retrieval]
      for j in range(0, len(p_values)):
        data2[k,i,j] = data1[j]
    
  for k in range(0,num_datasets):
    for i in range(0,len(S_values)):
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

#-----------------------------choose, only this section is to be touched------------------------------------------------------#

# put 1 for 0.7, 2 for 0.8, 3 for 0.9
retrieval = 1
frac_crit = 0.5
S_values = numpy.arange(1,7,1)
print(S_values)

#S=1
apf = [0.0, 0.4, 1.0, 0.4, 1.0, 0.05]
prob = [0.0, 0.05, 0.05, 0.2, 0.2, 0.05]
colors = ['black', 'green', 'blue', 'red' , 'purple', 'orange']
p_values = [numpy.arange(20,2000,10),numpy.arange(20,2000,10),numpy.arange(20,1020,10),numpy.arange(10,1010,10),numpy.arange(10,1010,10),numpy.arange(10,2010,10)]
num_datasets = [1,1,1,1,1,1]

S_values_ticks = [1,2,3,4,5,6]
#p_values_ticks = numpy.arange(-20,220,20)

#-------------------------------------------plot-------------------------------------------------------------------------------------------------#

fig = plt.figure(figsize=(5.5, 5))
ax= fig.add_subplot(1,1,1) 
#ax.set_title(r"$\displaystyle\ N=%s, Cm=%s, S=%s, U=%s, \beta=%s$"%(N, Cm, S, U, beta ), fontsize=13)

for i in range(numpy.size(apf)):
  p_c = compute_pc(S_values, p_values[i][:], num_datasets[i], apf[i], prob[i])
  if apf[i] == 0.05 and prob[i] == 0.05:
	  plt.errorbar(S_values, numpy.mean(p_c/Cm, axis = 0), yerr=numpy.std(p_c/Cm, axis= 0), marker= 'o', label=r"$\displaystyle\ a_{p}=%.2f, f=%.2f $"%(apf[i], prob[i]), color=colors[i])
  else:
	  plt.errorbar(S_values, numpy.mean(p_c/Cm, axis = 0), yerr=numpy.std(p_c/Cm, axis= 0), marker= 'o', color=colors[i])
#plt.axvline(5, ls='--', color='black')  
ax.set_xlabel('$S$') 
ax.set_ylabel(r'$\alpha_c$')
ax.set_xticks(S_values_ticks)
legend = ax.legend(loc='best')
ax.set_ylim(0,10)
#plt.savefig('p_c_corr_fS N=%s Cm=%s a=%.1f U=%s beta=%s.pdf'%(N, Cm, a, U, beta ), bbox_inches='tight')
#plt.savefig('p_c_corr_fS N=%s Cm=%s a=%.1f U=%s beta=%s.png'%(N, Cm, a, U, beta ), bbox_inches='tight')
plt.savefig("test")
plt.close()
