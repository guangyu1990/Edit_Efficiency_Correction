#!/usr/bin/python

import click
import re
import os
import math
import pandas as pd
import numpy as np
import scipy
from scipy import stats
from sklearn import linear_model

@click.command()
@click.option("--file","-f",help="input file")
@click.option("--effect","-e",help="sgRNA efficiency file")
@click.option("--growth","-g",help="growth file")
@click.option("--lib","-l",help="lib type")

def main(file, effect, growth, lib):
    Eratio = {}
    Growth = []
    ER = []
    tag = ['day5_day0_rep1', 'day5_day0_rep2', 'day10_day0_rep1', 'day10_day0_rep2', 'day10_day5_rep1', 'day10_day5_rep2', 'day5_day0', 'day10_day0', 'day10_day5', 'day5_lib', 'day10_lib']
    tidx = 1
    cutoff = 0.1
    if '_NGG' not in lib:
        cutoff = 0.2
    with open(effect) as ef:
        for eline in ef:
            einf = eline.strip().split('\t')
            if eline.startswith('sgRNA'):
                continue
            sg, eratio = einf[0], float(einf[tidx])
            if eratio < cutoff:
                continue
            Eratio[sg] = float(eratio)
            ER.append(float(eratio))
    er_mean = np.median(ER)
    with open(growth) as gf:
        for gline in gf:
            ginf = gline.strip().split('\t')
            if gline.startswith('tag'):
                continue
            grow = float(ginf[1].split('|')[1])
            Growth.append(float(grow))
    with open(file) as f:
        for line in f:
            inf = line.strip().split('\t')
            if line.startswith('sgRNA'):
                print line.strip()
                continue
            ge, dlib, d0r1, d0r2, d5r1, d5r2, d10r1, d10r2 = inf[0], inf[5], inf[6], inf[7], inf[8], inf[9], inf[10], inf[11]
            if float(d0r1)*float(d0r2)*float(d5r1)*float(d5r2)*float(d10r1)*float(d10r2) == 0:
                continue
            d5d0r1, d5d0r2, d10d0r1, d10d0r2, d10d5r1, d10d5r2, d5d0, d10d0, d10d5, d5lib, d10lib = float(inf[19]), float(inf[20]), float(inf[21]), float(inf[22]), float(inf[23]), float(inf[24]), float(inf[25]), float(inf[26]), float(inf[27]), float(inf[28]), float(inf[29])
            FC = [math.log(d5d0r1), math.log(d5d0r2), math.log(d10d0r1), math.log(d10d0r2), math.log(d10d5r1), math.log(d10d5r2), math.log(d5d0), math.log(d10d0), math.log(d10d5), math.log(d5lib), math.log(d10lib)]
            Score = []
            effect = 1
            if ge in Eratio and Eratio[ge] != 0:
                effect = Eratio[ge]
                ge = ge+'_'+str(Eratio[ge])
            #elif '_' not in ge:
            #    effect = er_mean
            #    ge = ge+'_'+str(er_mean)
            else:
                ge = ge+'_NA'   
           
            for i in range(0, len(FC)):
                s = (FC[i]-Growth[i])/effect
                #if '_NGG' not in lib:
                #    if s > 0:
                #        s = s/1.6
                #    if s < 0:
                #        s = s/1.6
                #if s > 2.5:
                #    s = 2.5
                #if s < -5:
                #    s = -5
                Score.append(str(s)) 
            sd5d0, sd10d0, sd10d5 = (float(Score[0])+float(Score[1]))/2, (float(Score[2])+float(Score[3]))/2, (float(Score[4])+float(Score[5]))/2  
            print ge+'\t'+'\t'.join(inf[1:19])+'\t'+'\t'.join(Score[0:6])+'\t'+str(sd5d0)+'\t'+str(sd10d0)+'\t'+str(sd10d5)+'\t'+'\t'.join(Score[9:len(Score)])

if __name__ == "__main__":
    main()
