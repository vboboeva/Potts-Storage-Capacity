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

# get constants from const_potts C++ file
lines = open( "const_potts.h", "r" ).readlines()

N = int( (lines[0].split())[2] )
S = int( (lines[1].split())[2] )
Cm = int( (lines[2].split())[2] )
p_fact = int( (lines[3].split())[2] )
Num_fact = int( (lines[4].split())[2] )
apf = float( (lines[5].split())[2] )
U = float( (lines[6].split())[2] )
beta = int( (lines[7].split())[2] )

print("N=",N)
print("Cm=",Cm)
print("S=",S)
print("U=",U)

#---------------------------------------function to compute p_c----------------------------------------------#

def compute_pc(a_values, p_values, num_datasets, S, apf):
  
  p_cc = numpy.zeros((num_datasets,len(a_values),len(p_values)))
  p_c = numpy.zeros((num_datasets,len(a_values)))
  i=0
  j=0
  data2 = numpy.zeros((num_datasets,len(a_values),len(p_values)))
  for k in range(0,num_datasets):
    data = (loadtxt('S%s_ap%.1f/storage_capacity_dataset_%s'%(S, apf, k)))
    for i in range(0,len(a_values)):
      data1 = data[(i)*len(p_values):(i+1)*len(p_values),retrieval]
      #print(numpy.size(data1))
      for j in range(0, len(p_values)):
        data2[k,i,j] = data1[j]
    
  for k in range(0,num_datasets):
    for i in range(0,len(a_values)):
      for j in range(0, len(p_values)):      
        if data2[k,i,j]>= 1:      
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
a_values = numpy.linspace(0.1,0.9,9)
print(a_values)

#S=2
apf = [0.0, 0.4, 1.0]
colors = ['green', 'blue', 'red' ]
p_values = [numpy.arange(20,100,1),numpy.arange(20,100,1),numpy.arange(20,100,1)]
num_datasets = [1,1,1]
a_values_ticks = numpy.linspace(0.0,1.00,11)
p_values_ticks = numpy.arange(-20,110,10)

#-------------------------------------------plot-------------------------------------------------------------------------------------------------#

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

fig = plt.figure(figsize=(5.5, 5))
ax= fig.add_subplot(1,1,1) 
#ax.set_title(r"$\displaystyle\ N=%s, Cm=%s, S=%s, U=%s, \beta=%s$"%(N, Cm, S, U, beta ), fontsize=13)

for i in range(numpy.size(apf)):
  print(apf[i])
  p_c = compute_pc(a_values, p_values[i][:], num_datasets[i], S, apf[i])
  plt.errorbar(a_values, numpy.mean(p_c, axis = 0), yerr=numpy.std(p_c, axis= 0), marker= 'o', label=r"$\displaystyle\ a_{p}=%.1f $"%(apf[i]), color=colors[i])

#data_ana =  (loadtxt("S2alphacCorrelated.txt"))
#plt.plot(data_ana[:,0], N*data_ana[:,1], linewidth = 3, linestyle = '--', color='black', label = r"analytic")
  
ax.set_xlabel('$a$', fontsize=18) 
ax.set_ylabel('$p_c$', fontsize=18)
ax.set_xticks(a_values_ticks)
#ax.set_yticks(p_values_ticks)
ax.tick_params(axis='x', labelsize=14)
ax.tick_params(axis='y', labelsize=14)
#plt.ylim(p_values_ticks[0], p_values_ticks[-1])
legend = ax.legend(loc='upper right', fontsize=16)
#plt.grid(True)
plt.savefig('p_c_storage N=%s Cm=%s S=%s U=%s beta=%s.png'%(N, Cm, S, U, beta ))
plt.close()
