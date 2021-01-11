#!/usr/bin/python

import click
import re
import os
import math
import pandas as pd
import numpy as np
import scipy
from scipy import stats
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import seaborn as sns

@click.command()
@click.option("--file","-f",help="input file in special format: sgRNA\tGene\tFunc\tClinic\tExon\tScore")

def main(file):
    Score, control, Var = {}, [], []
    Smean, Ssd = [], []
    with open(file) as f:
        for line in f:
            inf = line.strip().split('\t')
            if line.startswith('sgRNA'):
                print('sgRNA\tGene\tFunc\tClinic\tExon\tScore\tPvalue')
                continue
            sg, ge, func, scr_rep1, scr_rep2, score = inf[0], inf[1], inf[2], float(inf[21]), float(inf[22]), float(inf[26])
            if 'nontargeting' in ge:
                control.append(score)
            elif '_NA' in sg:
                continue
            elif 'BRCA' not in sg:
                continue
            scr_mean, scr_std = np.mean([scr_rep1, scr_rep2]), np.std([scr_rep1, scr_rep2])
            if scr_mean * scr_std != 0:
                Ssd.append(math.log(abs(scr_std*scr_std-scr_mean), 2)), Smean.append(math.log(abs(scr_mean), 2))
            Score[sg] = '\t'.join(inf[0:5])+'\t'+str(score)+'\t'+str(scr_mean)+'\t'+str(scr_std*scr_std)
            if 'tive' not in func:
                Var.append(score)
    mean, std, num = np.mean(control), np.std(control), len(control)
    vmean, vstd, vmedian = np.mean(Var), np.std(Var), np.median(Var)
    p = stats.ks_2samp(control, stats.norm.rvs(loc=mean, scale=std, size=num))[1]
    #print(str(mean)+'\t'+str(std)+'\t'+str(p))

    plt.figure(figsize=(5,3))        
    plt.rcParams['font.sans-serif'] = ['Arial']
    plt.rcParams['font.size'] = 8
    ax = plt.gca()
    z = np.polyfit(Smean,Ssd,1)
    a, b = float(z[0]), float(z[1])
    yf = []
    for x0 in np.sort(Smean):
        yf.append(a*(x0) + b)
    ax.plot(np.sort(Smean),yf,color="#999999",linestyle='-.')
    cor, pvalue = round(stats.pearsonr(Smean,Ssd)[0], 3), stats.pearsonr(Smean,Ssd)[1]
    ax.scatter(Smean,Ssd,color="#999999",s=5)
    ax.annotate('pearson_r='+str(cor),xy=(0,-10), horizontalalignment='left', verticalalignment='top',fontsize=8)
    ax.annotate('pearson_p='+str(pvalue),xy=(0,-12), horizontalalignment='left', verticalalignment='top',fontsize=8)
    plt.savefig(file+'.pdf',format='pdf',bbox_inches='tight')


    for s in Score:
        Score_inf = Score[s].split('\t')
        value = float(Score[s].split('\t')[5])
        value_mean = float(Score[s].split('\t')[6])
        value_var = float(Score[s].split('\t')[7])
        std_et = math.sqrt(abs(2**(a*math.log(abs(value_mean),2)+b)+value_mean))
        if value_var < value_mean:
            std_et = math.sqrt(abs(value_mean-2**(a*math.log(abs(value_mean),2)+b)))
        value = value/std_et
        pvalue = stats.gaussian_kde(stats.norm.rvs(loc=0, scale=1, size=10000)).integrate_box_1d(-np.Inf,value)
        print('\t'.join(Score_inf[0:6])+'\t'+str(pvalue))        



if __name__ == "__main__":
    main()
