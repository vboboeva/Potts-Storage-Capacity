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
  
  #print(len(zeta_values))
  #print(len(p_values))
  p_cc = numpy.zeros((num_datasets,len(zeta_values),len(p_values)))
  p_c = numpy.zeros((num_datasets,len(zeta_values)))
  i=0
  j=0
  data2 = numpy.zeros((num_datasets,len(zeta_values),len(p_values)))
  for k in range(0,num_datasets):
    data = (loadtxt('storage_capacity_dataset_%s'%(k)))
    for i in range(0,len(zeta_values)):
      data1 = data[(i)*len(p_values):(i+1)*len(p_values),retrieval]
      #print(data1)
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
  #print("p_c=", p_c) 
  return p_c
  
def compute_info(a_values, p_values, num_datasets, i):
	info, frac_retr, frac_retr_2ndpatt, frac_retr_par, mmean_retr, mmean_retr_2ndpatt, mmean_retr_par, sparsity = numpy.zeros((num_datasets, numpy.size(p_values))), numpy.zeros((num_datasets, numpy.size(p_values))), numpy.zeros((num_datasets, numpy.size(p_values))), numpy.zeros((num_datasets, numpy.size(p_values))), numpy.zeros((num_datasets, numpy.size(p_values))), numpy.zeros((num_datasets, numpy.size(p_values))), numpy.zeros((num_datasets, numpy.size(p_values))), numpy.zeros((num_datasets, numpy.size(p_values)))
	for k in range(0, num_datasets):
		data = (loadtxt('storage_capacity_dataset_%s'%(k)))
		pvals = data[(i)*len(p_values):(i+1)*len(p_values),0]
		info[k,:] = data[(i)*len(p_values):(i+1)*len(p_values),31]
		
		frac_retr[k,:] = data[(i)*len(p_values):(i+1)*len(p_values),retrieval]
		frac_retr_2ndpatt[k,:] = data[(i)*len(p_values):(i+1)*len(p_values),retrieval2]
		frac_retr_par[k,:] = data[(i)*len(p_values):(i+1)*len(p_values),retrieval3]
		
		mmean_retr[k,:] = data[(i)*len(p_values):(i+1)*len(p_values),1]
		mmean_retr_2ndpatt[k,:] = data[(i)*len(p_values):(i+1)*len(p_values),1+10]
		mmean_retr_par[k,:] = data[(i)*len(p_values):(i+1)*len(p_values),1+10+10]
		
		sparsity[k,:] = data[(i)*len(p_values):(i+1)*len(p_values),sp]
	
	#print(numpy.mean(info, axis=0))	
	return pvals, numpy.mean(info, axis=0), numpy.mean(frac_retr, axis=0), numpy.mean(frac_retr_2ndpatt, axis=0), numpy.mean(frac_retr_par, axis=0), numpy.mean(mmean_retr, axis=0), numpy.mean(mmean_retr_2ndpatt, axis=0), numpy.mean(mmean_retr_par, axis=0), numpy.mean(sparsity, axis=0)

#-----------------------------choose, only this section is to be touched------------------------------------------------------#

retrieval = 1+9
retrieval2 = retrieval+9
retrieval3 = retrieval+9+9
sp = 9
frac_crit = 0.5

color = 'green'
p_values = numpy.arange(10,2010,10)
zeta_values = numpy.array(( 0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1, 5, 10, 50, 100, 500))

plot_vals = numpy.array((0.001000, 0.01, 0.1, 1, 10))
num_datasets = 4

xsize=4*1.62
ysize=4

#----------------------------------------------------storage capacity

fig = plt.figure(figsize=(5.5, 5))
ax= fig.add_subplot(1,1,1) 
p_c = compute_pc(zeta_values, p_values, num_datasets)
plt.errorbar(zeta_values, numpy.mean(p_c/Cm, axis = 0), yerr=numpy.std(p_c/Cm, axis= 0), marker= 'o', color=color)

ax.set_xlabel(r'$\zeta$') 
ax.set_ylabel(r'$\alpha_c$')
plt.ylim(0, 8)
#plt.xlim(0,0.2)
#legend = ax.legend(loc='lower center')
#plt.grid(True)
plt.savefig('p_c_zeta retrieval=%s N=%s Cm=%s S=%s U=%s beta=%s.png'%(retrieval, N, Cm, S, U, beta ),bbox_inches='tight')
plt.savefig('p_c_zeta retrieval=%s N=%s Cm=%s S=%s U=%s beta=%s.pdf'%(retrieval, N, Cm, S, U, beta ),bbox_inches='tight')

#-------------------------------------------------information


k=0
meaninfovals = []
for j in range(0,len(zeta_values)):
	if (zeta_values[j] in plot_vals):
		fig = plt.figure(figsize=(xsize, ysize))
		ax = fig.add_subplot(1,1,1)
		pvals, info, frac_retr, frac_retr_2ndpatt, frac_retr_par, mmean_retr, mmean_retr_2ndpatt, mmean_retr_par, sparsity = compute_info(zeta_values, p_values, num_datasets, j)
		ax.plot(pvals/Cm, info/Cm, color='black', dashes=(3*k+1,2), label = '$\zeta=%s$'%zeta_values[j])
		#print(info/Cm)
		k=k+1
		
		ax.set_xlabel(r'$\alpha$') 
		ax.set_ylabel(r'$I/c_m$')
		entropy = (-(1-a)*numpy.log(1-a)+a*numpy.log(S/a))/(numpy.log(2)*Cm)
		#print(entropy)
		ax.axhline(entropy, lw=2.5, label = '$entropy$')
		plt.ylim(0, 0.0036)
		legend = ax.legend(loc='best')
		plt.legend()
		plt.savefig('info N=%s Cm=%s S=%s U=%s beta=%s zeta=%3f.png'%(N, Cm, S, U, beta, zeta_values[j]),bbox_inches='tight')
		plt.savefig('info N=%s Cm=%s S=%s U=%s beta=%s zeta=%3f.pdf'%(N, Cm, S, U, beta, zeta_values[j]),bbox_inches='tight')


	pvals, info, frac_retr, frac_retr_2ndpatt, frac_retr_par, mmean_retr, mmean_retr_2ndpatt, mmean_retr_par, sparsity = compute_info(zeta_values, p_values, num_datasets, j)
	x = numpy.where(frac_retr<0.1)
	#print(x)
	meaninfovals = numpy.append(meaninfovals,  numpy.mean(info[x]))

# plot the residual
fig = plt.figure(figsize=(xsize, ysize))
ax = fig.add_subplot(1,1,1)	
#print(meaninfovals/Cm)
ax.plot(zeta_values, meaninfovals/Cm, 'ko-')
ax.set_xlabel(r'$\zeta$') 
ax.set_ylabel(r'$I_r/c_m$')
k=0
for j in range(0,len(zeta_values)):
	if (zeta_values[j] in plot_vals):
		plt.axvline(zeta_values[j], dashes=(3*k+1,2), color='black')
		k=k+1
plt.xscale('log')
plt.savefig('res_info_fzeta N=%s Cm=%s S=%s U=%s beta=%s.png'%(N, Cm, S, U, beta ),bbox_inches='tight')
plt.savefig('res_info_fzeta N=%s Cm=%s S=%s U=%s beta=%s.pdf'%(N, Cm, S, U, beta ),bbox_inches='tight')	


#-----------------------------------------------------------------------fraction retrieved

k=0
for j in range(0,len(zeta_values)):
	if (zeta_values[j] in plot_vals):
		fig = plt.figure(figsize=(xsize, ysize))
		ax2 = fig.add_subplot(1,1,1)
		ax3 = ax2.twinx()

		pvals, info, frac_retr, frac_retr_2ndpatt, frac_retr_par,mmean_retr, mmean_retr_2ndpatt, mmean_retr_par, sparsity = compute_info(zeta_values, p_values, num_datasets, j)
		ax2.plot(pvals/Cm, frac_retr, color='green', dashes=(3*k+1,2), label = '$\zeta=%s$'%zeta_values[j])
		ax3.plot(pvals/Cm, frac_retr_2ndpatt, color='blue', dashes=(3*k+1,2))
		ax3.plot(pvals/Cm, frac_retr_par, color='red', dashes=(3*k+1,2))
		#ax3.plot(pvals/Cm, sparsity, ':', color=colors[k])
		k=k+1
	
		ax2.set_xlabel(r'$\alpha$') 
		ax2.set_ylabel(r'\bf{fraction retrieved}', color='green')
		ax3.set_ylabel(r'\bf{correlated pattern (blue)/ parent (red)}')
		ax2.set_ylim(0,1)
		ax3.set_ylim(0,1)
		ax2.tick_params(axis='x')
		ax2.tick_params(axis='y', colors='green')
		#ax3.tick_params(axis='y', colors='red')
		
		#plt.xlim(0, 4)
		legend = ax2.legend(loc='center right')
		plt.savefig('frac retr retrieval=%s N=%s Cm=%s S=%s U=%s beta=%s zeta=%3f.png'%(retrieval, N, Cm, S, U, beta, zeta_values[j] ),bbox_inches='tight')
		plt.savefig('frac retr retrieval=%s N=%s Cm=%s S=%s U=%s beta=%s zeta=%3f.pdf'%(retrieval, N, Cm, S, U, beta, zeta_values[j] ),bbox_inches='tight')
	

#-----------------------------------------------------------------------mean overlap

k=0
for j in range(0,len(zeta_values)):
	if (zeta_values[j] in plot_vals):
		fig = plt.figure(figsize=(xsize, ysize))
		ax2 = fig.add_subplot(1,1,1)
		ax3 = ax2.twinx()
		
		pvals, info, frac_retr, frac_retr_2ndpatt, frac_retr_par, mmean_retr, mmean_retr_2ndpatt, mmean_retr_par, sparsity = compute_info(zeta_values, p_values, num_datasets, j)
		ax2.plot(pvals/Cm, mmean_retr, color='green', dashes=(3*k+1,2), label = '$\zeta=%s$'%zeta_values[j])
		ax3.plot(pvals/Cm, mmean_retr_2ndpatt, color='blue', dashes=(3*k+1,2))
		ax3.plot(pvals/Cm, mmean_retr_par, color='red', dashes=(3*k+1,2))
		k=k+1
	
		ax2.set_xlabel(r'$\alpha$') 
		ax2.set_ylabel(r'\bf{mean overlap}', color='green')
		ax3.set_ylabel(r'\bf{correlated pattern (blue)/ parent (red)}')
		ax2.set_ylim(0,1)
		ax3.set_ylim(0,1)
		ax2.tick_params(axis='x')
		ax2.tick_params(axis='y', colors='green')
		#ax3.tick_params(axis='y', colors='red')

		#plt.xlim(0, 4)
		legend = ax2.legend(loc='center right')
		plt.savefig('mean retrieval=%s N=%s Cm=%s S=%s U=%s beta=%s zeta=%3f.png'%(retrieval, N, Cm, S, U, beta, zeta_values[j] ),bbox_inches='tight')
		plt.savefig('mean retrieval=%s N=%s Cm=%s S=%s U=%s beta=%s zeta=%3f.pdf'%(retrieval, N, Cm, S, U, beta, zeta_values[j] ),bbox_inches='tight')
