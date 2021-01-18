#!/usr/bin/python

import click
import re
import os
import math
import pandas as pd
import numpy as np
import scipy
from scipy import stats
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, auc
#from sklearn import cross_validation

@click.command()
@click.option("--file","-f",help="data file")

def main(file):
    Y_true, Y_score_our, Y_score_mageck = [], [], []
    #EF = [0, 0.001, 0.005, 0.01, 0.05, 0.1]
    EF = [0.1, 0.2]
    for e in EF:
        y_true, y_score_our, y_score_mageck = [], [], []
        Sample = [] 
        with open(file) as f:
            for line in f:
                if line.startswith('sgRNA'):
                    continue
                inf = line.strip().split('\t')
                if inf[5] == 'NA':
                    continue
                sg, tag, fdr, mbeta, mpvalue, ef, score, es, rra = inf[0], inf[9], float(inf[5]), float(inf[7]), float(inf[8]), inf[10], float(inf[4]), float(inf[11]), float(inf[13])
                tag_value = 0
                if tag == 'LOF':
                    tag_value = 1
                if ef == 'NA':
                    ef = 0
                if float(ef) == 0:
                    continue
                if float(es) < e:
                    continue     
                y_true.append(tag_value)
                Sample.append(sg)
                y_score_our.append(fdr)
                y_score_mageck.append(1-mpvalue)

        Y_true.append(y_true)
        Y_score_our.append(y_score_our)
        Y_score_mageck.append(y_score_mageck)
    FPR, TPR, AUC, SPC = {}, {}, {}, {}
    for i in range(0, len(EF)):            
        fpr_our, tpr_our, thresholds_our = roc_curve_self(Y_true[i], Y_score_our[i])
        spc_our = get_specificity(fpr_our)
        auc_our = auc(fpr_our, tpr_our)
        fpr_mageck, tpr_mageck, thresholds_mageck = roc_curve_self(Y_true[i], Y_score_mageck[i])
        spc_mageck = get_specificity(fpr_mageck)
        auc_mageck = auc(fpr_mageck, tpr_mageck)
        print '#raw efficiency cut off:'+str(EF[i])
        print 'AUC of our\t'+str(auc_our)
        print 'AUC of mageck\t'+str(auc_mageck)
        FPR.setdefault(EF[i], [])
        FPR[EF[i]] = [fpr_our, fpr_mageck]
        TPR.setdefault(EF[i], [])
        TPR[EF[i]] = [tpr_our, tpr_mageck]
        AUC.setdefault(EF[i], [])
        AUC[EF[i]] = [auc_our, auc_mageck]
        SPC.setdefault(EF[i], [])
        SPC[EF[i]] = [spc_our, spc_mageck]

    draw_roccurve(SPC, TPR, AUC, EF)

def roc_curve_self(ytrue, yscore):
    score_tmp = []
    for yn in yscore:
        score_tmp.append(yn)
    yscore.sort()
    true_tmp = ytrue
    ytrue = []
    for yi in range(0, len(yscore)):
        for yj in range(0, len(score_tmp)):
            if yscore[yi] == score_tmp[yj]:
                ytrue.append(true_tmp[yj])  
                break
    #print len(ytrue)

    fpr, tpr, thred = [], [], []
    for scr in yscore:
        ypred = []
        tp, fp, tn, fn = 0, 0, 0, 0
        for sc in yscore:
            ytg = 0
            if sc >= scr:
                ytg = 1     
            ypred.append(ytg)
        for n in range(0,len(ytrue)):
            if ytrue[n]*ypred[n] == 1:
                tp += 1
            else:
                if ytrue[n] == 1:
                    fn += 1
                elif ypred[n] == 1:
                    fp += 1
                else:
                    tn += 1
        fpr_n = float(fp)/(fp+tn)
        tpr_n = float(tp)/(tp+fn)
        fpr.append(fpr_n)
        tpr.append(tpr_n)
        thred.append(scr)
    return fpr, tpr, thred
        

def get_specificity(fpr):
    spc = []
    for f in range(0, len(fpr)):
        spc.append(1-fpr[f])
    return spc

def draw_roccurve(SPC, TPR, AUC, EF):
    plt.figure(figsize=(5,5))
    plt.rcParams['font.sans-serif'] = ['Arial']
    plt.rcParams['font.size'] = 8
    colors = ["#2166AC", "#4393C3", "#92C5DE", "#D1E5F0", "#FDDBC7", "#F4A582", "#D6604D", "#B2182B"]
    colors = [ "#D1E5F0", "#92C5DE", "#4393C3", "#2166AC", "#FDDBC7", "#F4A582", "darkorange", "#D55E00"]
    colors = [ "#92C5DE", "#4393C3", "#2166AC", "#F4A582", "darkorange", "#D55E00"]
    colors = [ "#4393C3", "#2166AC", "darkorange", "#D55E00"]
    lw = 1
    
    for j in range(0, len(EF)):
        plt.plot(SPC[EF[j]][0], TPR[EF[j]][0], color=colors[j],lw=lw, label='Efficiency-corrected '+str(EF[j])+' AUC = %0.2f' % AUC[EF[j]][0])
        plt.plot(SPC[EF[j]][1], TPR[EF[j]][1], color=colors[j+2],lw=lw, label='MAGeCK '+str(EF[j])+' AUC = %0.2f' % AUC[EF[j]][1])
    plt.plot([0, 1], [1, 0], color='grey', lw=lw, linestyle='--')
    plt.xlim([0, 1])
    plt.ylim([0, 1])
    plt.xlabel('Specificity')
    plt.ylabel('Sensitivity')
    plt.legend(loc="lower left", frameon=False)
    #plt.title('EfficiencyCorrectedModel ROC')
    ax=plt.gca()
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    plt.savefig('compare_roc.pdf',format='pdf',bbox_inches='tight')


if __name__ == "__main__":
    main()
