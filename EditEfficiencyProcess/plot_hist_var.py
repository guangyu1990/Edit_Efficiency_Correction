#!/usr/bin/env python
import os
import sys
import numpy as np
from mpl_toolkits.axes_grid1 import host_subplot
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
import click
import seaborn as sns
import math

file = sys.argv[1]
outfile = sys.argv[2]+'.png'
x = []
xrb, xrp, xr, xw, xt, xok = [], [], [], [], [], []
#yrb, yrp, yr, yw, yt, yok = [], [], [], [], [], []
crb, crp, cr, cw, ct, cok = [], [], [], [], [], []
ff = open(file)
j = 1
for line in ff:
    if line.startswith('#'):
        continue
    items = line.rstrip().split('\t')
    vtype = items[7]
    value = float(items[6])
    x.append(value)
    if vtype == 'target':
        xt.append(value)
        xok.append(value)
    elif vtype == 'window':
        xw.append(value)
        xok.append(value)
    elif vtype == 'wrong pos':
        xrp.append(value)
        xr.append(value)
    elif vtype == 'wrong base':
        xrb.append(value)
        xr.append(value)    

ff.close()

plt.figure(figsize=(12,4))
plt.rcParams['font.sans-serif'] = ['Arial']
plt.rcParams['font.size'] = 12

#ax = plt.gca()
ax = plt.subplot(111)
hand, lab = [], []
n_bins = 500


#b1=ax.hist(xr, n_bins, histtype='bar', color="#999999", edgecolor = 'k', label="all wrong")
#b2=ax.hist(xrb, n_bins, histtype='bar', color="#0072B2", edgecolor = 'k', label="wrong base")
#b3=ax.hist(xrp, n_bins, histtype='bar', color="#56B4E9", edgecolor = 'k', label="wrong pos")
#b4=ax.hist(xok, n_bins, histtype='bar', color="red", edgecolor = 'k', label="all ok")
#b5=ax.hist(xw, n_bins, histtype='bar', color="#E69F00", edgecolor = 'k', label="window")
b6=ax.hist(xt, n_bins, histtype='bar', color="#D55E00", edgecolor = 'k', label="target")

'''
sns.set_style("white")
for i in range(0,len(X)):
    if len(X[i]) > 1 :
        #B[i] = sns.distplot(X[i], color=C[i], label=F[i]).legend()
        B[i] = sns.kdeplot(X[i], shade=True, color=C[i], label=F[i]).legend()
'''
#ax.set_xlim([0, 1000])
ax.legend()
plt.tight_layout()
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.set_xlabel("freq")
ax.set_ylabel("num of var")
plt.savefig(outfile,format='png',bbox_inches='tight')
