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
@click.option("--file","-f",help="input file")
@click.option("--outfile","-o",help="output file")

def main(file, outfile):
    Grow = {}
    tag = ['day5_day0_rep1', 'day5_day0_rep2', 'day10_day0_rep1', 'day10_day0_rep2', 'day10_day5_rep1', 'day10_day5_rep2', 'day5_day0', 'day10_day0', 'day10_day5', 'day5_lib', 'day10_lib']
    colors = ["#CC79A7", "#9999CC", "#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00","#CC6666", "#66CC99"]
    with open(file) as f:
        for line in f:
            inf = line.strip().split('\t')
            if line.startswith('sgRNA'):
                #print line.strip()
                for i in range(19,30):
                    Grow.setdefault(inf[i], [])
                continue
            lib, d0r1, d0r2, d5r1, d5r2, d10r1, d10r2 = inf[5], inf[6], inf[7], inf[8], inf[9], inf[10], inf[11]
            if float(d0r1)*float(d0r2)*float(d5r1)*float(d5r2)*float(d10r1)*float(d10r2) == 0:
                continue
            d5d0r1, d5d0r2, d10d0r1, d10d0r2, d10d5r1, d10d5r2, d5d0, d10d0, d10d5, d5lib, d10lib = float(inf[19]), float(inf[20]), float(inf[21]), float(inf[22]), float(inf[23]), float(inf[24]), float(inf[25]), float(inf[26]), float(inf[27]), float(inf[28]), float(inf[29])
            value = [math.log(d5d0r1), math.log(d5d0r2), math.log(d10d0r1), math.log(d10d0r2), math.log(d10d5r1), math.log(d10d5r2), math.log(d5d0), math.log(d10d0), math.log(d10d5), math.log(d5lib), math.log(d10lib)]
            for j in range(0,len(tag)):
                Grow[tag[j]].append(value[j])
            #print '\t'.join(inf[0:19])+'\t'+str(d5d0r1)+'\t'+str(d5d0r2)+'\t'+str(d10d0r1)+'\t'+str(d10d0r2)+'\t'+str(d10d5r1)+'\t'+str(d10d5r2)+'\t'+str(d5d0)+'\t'+str(d10d0)+'\t'+str(d10d5)+'\t'+str(d5lib)+'\t'+str(d10lib)
    print 'tag\tgrowth_rate\tstd\tp-value'
    tmp = [2,3,7]
    for m in range(0,len(tag)):
        array = np.array(Grow[tag[m]])
        mean, std, median = np.mean(array), np.std(array), np.median(array)
        ax = plt.gca()
        sns.set(color_codes=True)
        sns.set_style("white")
        sns.kdeplot(array, shade=True, color=colors[m], label=tag[m]).legend()
        array_std = (array-np.mean(array))/np.std(array)
        nm_out = stats.ks_2samp(array, stats.norm.rvs(loc=mean, scale=std, size=len(array)))
        #nm_out = scipy.stats.kstest(array_std, 'norm')[1]
        #sns.kdeplot(array_std, shade=True, color=colors[m], label=tag[m]).legend()
        ax.legend()
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.set_xlabel("ln(growth rate)")
        ax.set_ylabel("num of SNV")
        plt.savefig(outfile+'.png',format='png',bbox_inches='tight')
        array = (array-np.mean(array))/np.std(array)
        nm_out = scipy.stats.kstest(array, 'norm')[1]
        print tag[m]+'\t'+str(mean)+'|'+str(median)+'\t'+str(std)+'\t'+str(nm_out)

def func(x,a,u,sig):
    return a*np.exp(-(x - u) ** 2 / (2 * sig ** 2)) / (sig * math.sqrt(2 * math.pi))


if __name__ == "__main__":
    main()
