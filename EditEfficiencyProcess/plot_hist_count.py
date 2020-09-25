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
ff = open(file)
j = 1
for line in ff:
    if '_' in line:
        items = line.rstrip().split('\t')
        value = float(items[1])
        x.append(value)
ff.close()

plt.figure(figsize=(12,4))
plt.rcParams['font.sans-serif'] = ['Arial']
plt.rcParams['font.size'] = 12

ax = plt.gca()
hand, lab = [], []
n_bins = 5000

b=ax.hist(x, n_bins, histtype='bar', color="#999999")
'''
sns.set_style("white")
for i in range(0,len(X)):
    if len(X[i]) > 1 :
        #B[i] = sns.distplot(X[i], color=C[i], label=F[i]).legend()
        B[i] = sns.kdeplot(X[i], shade=True, color=C[i], label=F[i]).legend()
'''
ax.set_xlim([0, 1000])
ax.legend()
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.set_xlabel("reads")
ax.set_ylabel("num of SNV")
plt.savefig(outfile,format='png',bbox_inches='tight')
