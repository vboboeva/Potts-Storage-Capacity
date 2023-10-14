# -*- coding: utf-8 -*-
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
from matplotlib import rc
from scipy import stats
from itertools import combinations
from ast import literal_eval
from matplotlib.colors import ListedColormap, NoNorm
from pylab import rcParams


# the axes attributes need to be set before the call to subplot
plt.rc('text', usetex=True)
plt.rc('font', family='serif', weight='bold', size=18)
plt.rc('xtick.major', size=5, pad=7)
plt.rc('ytick.major', size=5, pad=7)
rcParams['text.latex.preamble'] = [r'\usepackage{sfmath} \boldmath']
plt.rc('xtick', labelsize=18)
plt.rc('ytick', labelsize=18)


lines = open( "const_potts.h", "r" ).readlines()

N = int( (lines[0].split())[2] )
S = int( (lines[1].split())[2] )
p = int( (lines[2].split())[2] )
a = float( (lines[3].split())[2] )
apf = 1.0#float( (lines[6].split())[2] )
zeta = 0.000
N_fact=1000

#for low values of f, many children are globally null, so I discard them. The non-null number (pp) fluctuates.
#if f=0.03 pp=990
#if f=0.01 pp=784
#if f=0.05 pp=1000
pp=1000

for Num_fact in [150]:
	for p_fact in [30]:
		
		print(p_fact, Num_fact)		
		fig = plt.figure(figsize=(5, 5))

		Cmn = (loadtxt("Cmn1_file_a%.2f_apf%.2f_pfact%d_p%d_Nfact%d_Numfact%d_zeta%.3f"%(a, apf, p_fact, p, N_fact, Num_fact, zeta)))[numpy.triu_indices(pp, 1)]
		print(numpy.shape(Cmn))
		b = Cmn/a#[numpy.nonzero(Cmn)]/a
		print(numpy.mean(b))
		#plt.hist(b, bins='fd', histtype='step', normed=1, color='green', alpha=1, zorder=1, label = r'$C_{\mu \nu}$')
		hist, bin_edges = numpy.histogram(b, bins='sqrt', density=True)
		bin_centers = numpy.asarray([(bin_edges[i]+bin_edges[i+1])/2 for i in range(numpy.size(bin_edges)-1)])
		indices = numpy.where(hist>0)
		hist_new=hist[indices]
		bins_new=bin_centers[indices]
		plt.fill_between(bins_new,0,hist_new, alpha=0.5,color='green',label = r'$C_{\mu \nu}$')
		
		
		Cij = (loadtxt("Cij1_file_a%.2f_apf%.2f_pfact%d_p%d_Nfact%d_Numfact%d_zeta%.3f"%(a, apf, p_fact, p, N_fact, Num_fact, zeta)))[numpy.triu_indices(N, 1)]
		#print(numpy.nonzero(Cij))
		b = Cij/a#[numpy.nonzero(Cij)]/a
		print(numpy.mean(b))
		#plt.hist(b, bins='auto', histtype='step', normed=1, color='red', alpha=1, zorder=2, label = r'$C_{ij}$')
		hist, bin_edges = numpy.histogram(b, bins='sqrt', density=True)
		bin_centers = numpy.asarray([(bin_edges[i]+bin_edges[i+1])/2 for i in range(numpy.size(bin_edges)-1)])
		indices = numpy.where(hist>0)
		hist_new=hist[indices]
		bins_new=bin_centers[indices]
		plt.fill_between(bins_new,0,hist_new, alpha=0.5,color='red', label = r'$C_{ij}$')
		
		#plt.xticks(fontsize=14)
		#plt.yticks(fontsize=14)
		plt.ylim(0.01,4000)
		plt.xlim(0,0.4)
		plt.axvline(a/S, color='black')
		plt.yscale('log', nonposy='clip')
		#plt.xlabel('correlation', fontsize=18)
		if (zeta == 0.001): 
			#plt.legend()
			fig.suptitle('$f=%.2f$'%(float(p_fact)/float(p)))
		if (p_fact == 50):
			plt.annotate("$\zeta=%s$"%zeta, (0.25,1))
		fig.savefig("distribs_a%.1f_apf%.1f_pfact%d_p%d_Nfact%d_Numfact%d_zeta%.3f.pdf"%(a, apf, p_fact, p, N_fact, Num_fact, zeta), bbox_inches='tight')
		fig.savefig("distribs_a%.1f_apf%.1f_pfact%d_p%d_Nfact%d_Numfact%d_zeta%.3f.png"%(a, apf, p_fact, p, N_fact, Num_fact, zeta), bbox_inches='tight')
