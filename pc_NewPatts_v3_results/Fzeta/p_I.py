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
#---------------------------------------function to compute p_c----------------------------------------------#

def compute_pc(zeta_values, p_values, num_datasets):
  
  print(len(zeta_values))
  print(len(p_values))
  p_cc = numpy.zeros((num_datasets,len(zeta_values),len(p_values)))
  p_c = numpy.zeros((num_datasets,len(zeta_values)))
  i=0
  j=0
  data2 = numpy.zeros((num_datasets,len(zeta_values),len(p_values)))
  for k in range(0,num_datasets):
    data = (loadtxt('apf%.1f_f%.2f/storage_capacity_dataset_%s'%(apf, f, k)))
    for i in range(0,len(zeta_values)):
      data1 = data[(i)*len(p_values):(i+1)*len(p_values),retrieval]
      #print(numpy.size(data1))
      for j in range(0, len(p_values)):
        data2[k,i,j] = data1[j]
    
  for k in range(0,num_datasets):
    for i in range(0,len(zeta_values)):
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
		info = data[(i)*len(p_values):(i+1)*len(p_values),7]
		frac_retr = data[(i)*len(p_values):(i+1)*len(p_values),retrieval]
		frac_retr_2ndpatt = data[(i)*len(p_values):(i+1)*len(p_values),retrieval2]
		sparsity = data[(i)*len(p_values):(i+1)*len(p_values),sp]
	return pvals, info, frac_retr, frac_retr_2ndpatt, sparsity

#-----------------------------choose, only this section is to be touched------------------------------------------------------#

# put 1 for 0.7, 2 for 0.8, 3 for 0.9, 4 for 0.7corr, 5 for 0.8corr, 6 for 0.9corr, 7 for info, 9 for sparsity
retrieval = 2
retrieval2 = retrieval+3
sp = 9
frac_crit = 0.5

color = 'green'
p_values = numpy.arange(10,2010,10)
zeta_values = numpy.array(( 0.001, 0.002, 0.005, 0.008, 0.009, 0.01, 0.011, 0.012, 0.015, 0.02, 0.05, 0.1, 0.15, 0.2, 0.5, 1, 2, 3, 4, 5, 7, 50, 100, 500))
#zeta_values = numpy.array(( 0.001000, 0.002000, 0.005000, 0.008000, 0.009000, 0.010000, 0.011000, 0.012000, 0.015000, 0.020000, 0.050000))

plot_vals = numpy.array((0.001000, 0.01, 0.1, 1, 10, 100))
num_datasets = 1

xsize=4*1.62
ysize=4

#----------------------------------------------------storage capacity

fig = plt.figure(figsize=(5.5, 5))
ax= fig.add_subplot(1,1,1) 
p_c = compute_pc(zeta_values, p_values, num_datasets)
plt.errorbar(zeta_values, numpy.mean(p_c/Cm, axis = 0), yerr=numpy.std(p_c/Cm, axis= 0), marker= 'o', color=color)
plt.xscale('log')
ax.set_xlabel(r'$\zeta$') 
ax.set_ylabel(r'$\alpha_c$')
plt.ylim(0, 8)
#plt.xlim(0,0.2)
#legend = ax.legend(loc='lower center')
#plt.grid(True)
plt.savefig('p_c_zeta retrieval=%s N=%s Cm=%s S=%s U=%s beta=%s.png'%(retrieval, N, Cm, S, U, beta ),bbox_inches='tight')
plt.savefig('p_c_zeta retrieval=%s N=%s Cm=%s S=%s U=%s beta=%s.pdf'%(retrieval, N, Cm, S, U, beta ),bbox_inches='tight')

#-------------------------------------------------information

fig = plt.figure(figsize=(xsize, ysize))
ax = fig.add_subplot(1,1,1)
k=0
meaninfovals = []
for j in range(0,len(zeta_values)):
	if (zeta_values[j] in plot_vals):
		pvals, info, frac_retr, frac_retr_2ndpatt, sparsity = compute_info(zeta_values, p_values, num_datasets, j)
		ax.plot(pvals/Cm, info/Cm, color='black', dashes=(3*k+1,2), label = '$\zeta=%s$'%zeta_values[j])
		print(info/Cm)
		k=k+1
	pvals, info, frac_retr, frac_retr_2ndpatt, sparsity = compute_info(zeta_values, p_values, num_datasets, j)
	x = numpy.where(frac_retr<0.1)
	#print(x)
	meaninfovals = numpy.append(meaninfovals,  numpy.mean(info[x]))
ax.set_xlabel(r'$\alpha$') 
ax.set_ylabel(r'$I/c_m$')
entropy = (-(1-a)*numpy.log(1-a)+a*numpy.log(S/a))/(numpy.log(2)*Cm)
print(entropy)
ax.axhline(entropy, lw=2.5, label = '$entropy$')
plt.ylim(0, 0.0036)
legend = ax.legend(loc='best')
plt.legend()
plt.savefig('info N=%s Cm=%s S=%s U=%s beta=%s.png'%(N, Cm, S, U, beta ),bbox_inches='tight')
plt.savefig('info N=%s Cm=%s S=%s U=%s beta=%s.pdf'%(N, Cm, S, U, beta ),bbox_inches='tight')

# plot the residual
fig = plt.figure(figsize=(xsize, ysize))
ax = fig.add_subplot(1,1,1)	
print(meaninfovals/Cm)
ax.plot(zeta_values, meaninfovals/Cm, 'ko-')
ax.set_xlabel(r'$\zeta$') 
ax.set_ylabel(r'$I_r/c_m$')
#k=0
#for j in range(0,len(zeta_values)):
#	if (zeta_values[j] in plot_vals):
#		plt.axvline(zeta_values[j], dashes=(3*k+1,2), color='black')
#		k=k+1
plt.xscale('log')
plt.savefig('res_info_fzeta N=%s Cm=%s S=%s U=%s beta=%s.png'%(N, Cm, S, U, beta ),bbox_inches='tight')
plt.savefig('res_info_fzeta N=%s Cm=%s S=%s U=%s beta=%s.pdf'%(N, Cm, S, U, beta ),bbox_inches='tight')	


#-----------------------------------------------------------------------fraction retrieved

fig = plt.figure(figsize=(xsize, ysize))
ax2 = fig.add_subplot(1,1,1)
ax3 = ax2.twinx()
k=0
for j in range(0,len(zeta_values)):
	if (zeta_values[j] in plot_vals):
		pvals, info, frac_retr, frac_retr_2ndpatt, sparsity = compute_info(zeta_values, p_values, num_datasets, j)
		ax2.plot(pvals/Cm, frac_retr, color='green', dashes=(3*k+1,2), label = '$\zeta=%s$'%zeta_values[j])
		ax3.plot(pvals/Cm, frac_retr_2ndpatt, color='red', dashes=(3*k+1,2))
		#ax3.plot(pvals/Cm, sparsity, ':', color=colors[k])
		k=k+1
	
ax2.set_xlabel(r'$\alpha$') 
ax2.set_ylabel(r'\bf{fraction retrieved}', color='green')
ax3.set_ylabel(r'\bf{fraction retrieved another pattern}', color='red')
ax2.set_ylim(0,1)
ax3.set_ylim(0,1)
ax2.tick_params(axis='x')
ax2.tick_params(axis='y', colors='green')
ax3.tick_params(axis='y', colors='red')

#plt.xlim(0, 4)
legend = ax2.legend(loc='center right')
plt.legend()
plt.savefig('frac retr retrieval=%s N=%s Cm=%s S=%s U=%s beta=%s.png'%(retrieval, N, Cm, S, U, beta ),bbox_inches='tight')
plt.savefig('frac retr retrieval=%s N=%s Cm=%s S=%s U=%s beta=%s.pdf'%(retrieval, N, Cm, S, U, beta ),bbox_inches='tight')

#----------------------------------------------------------------------TOTAL retrieved
#color='black'
#fig = plt.figure(figsize=(xsize, ysize))
#ax2 = fig.add_subplot(1,1,1)
#k=0
#for j in range(0,len(zeta_values)):
	#if (zeta_values[j] in plot_vals):
		#pvals, info, frac_retr, frac_retr_2ndpatt, sparsity = compute_info(zeta_values, p_values, num_datasets, j)
		#ax2.plot(pvals/Cm, frac_retr+frac_retr_2ndpatt, color=color, dashes=(3*k+1,2))
		#k=k+1
	
#ax2.set_xlabel(r'$\alpha$') 
#ax2.set_ylabel(r'\bf{total retrieved}')
#ax2.set_ylim(0, 1)
#legend = ax2.legend(loc='center right')
#plt.savefig('total retr retrieval=%s N=%s Cm=%s S=%s U=%s beta=%s.png'%(retrieval, N, Cm, S, U, beta ),bbox_inches='tight')
#plt.savefig('total retr retrieval=%s N=%s Cm=%s S=%s U=%s beta=%s.pdf'%(retrieval, N, Cm, S, U, beta ),bbox_inches='tight')

#---------------------------------------------------------------------sparsity

color='black'
fig = plt.figure(figsize=(xsize, ysize))
ax2 = fig.add_subplot(1,1,1)
k=0
for j in range(0,len(zeta_values)):
	if (zeta_values[j] in plot_vals):
		pvals, info, frac_retr, frac_retr_2ndpatt, sparsity = compute_info(zeta_values, p_values, num_datasets, j)
		ax2.plot(pvals/Cm, sparsity, color=color, dashes=(3*k+1,2), label = '$\zeta=%s$'%zeta_values[j])
		k=k+1
	
ax2.set_xlabel(r'$\alpha$') 
ax2.set_ylabel(r'\bf{sparsity}')
plt.ylim(0, 1)
#legend = ax2.legend(loc='center right', fontsize=16)
plt.savefig('sparsity N=%s Cm=%s S=%s U=%s beta=%s.png'%(N, Cm, S, U, beta ),bbox_inches='tight')
plt.savefig('sparsity N=%s Cm=%s S=%s U=%s beta=%s.pdf'%(N, Cm, S, U, beta ),bbox_inches='tight')

