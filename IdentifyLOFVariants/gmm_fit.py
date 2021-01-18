#!/usr/bin/python

import click
import re
import os
import math
import pandas as pd
import numpy as np
import scipy
from scipy import stats
from sklearn import mixture
import matplotlib.pyplot as plt
import seaborn as sns

@click.command()
@click.option("--file","-f",help="input file")
@click.option("--num","-n",type=int,help="the number of components")
@click.option("--output","-o",help="output")

def main(file, num, output):
    data = pd.read_table(file, sep='\t', header='infer')
    data = data[data["day10_day0"] <= 3]
    score = data["day10_day0"].values
    #print data
    samp = len(score)
    SG, GE, FUNC, CLINC, EXON = data["sgRNA"].tolist(), data["Gene"].tolist(), data["Func"].tolist(), data["Clinic"].tolist(), data["Exon"].tolist()
    outfile = output+'_pred.txt'
    outfig = output+'.pdf'

    # Normal fit (for benchmark)
    score_grid = np.arange(min(score),max(score)+0.001,0.001);
    mean = np.mean(score)
    sigma2 = np.var(score)
    normal = stats.norm.pdf(score_grid, mean, math.sqrt(sigma2));

    # GMM fit
    gmm = mixture.GaussianMixture(n_components=num)
    gmm.fit(np.reshape(score,(samp,1)))
    aic, bic = gmm.aic(np.reshape(score,(samp,1))), gmm.bic(np.reshape(score,(samp,1)))
    
    #gauss_mixt = np.exp(gmm.score(np.reshape(score_grid,(len(score_grid),1))))
    gauss_mixt = np.array([p * stats.norm.pdf(score_grid, mu, sd) for mu, sd, p in zip(gmm.means_.flatten(), np.sqrt(gmm.covariances_.flatten()), gmm.weights_)])
    #gauss_mixt = np.sum(gauss_mixt, axis = 0)
    
    # Output
    gmm_mean = 0
    print 'weight\tmean\tstd'
    for mu, sd, p in zip(gmm.means_.flatten(), np.sqrt(gmm.covariances_.flatten()), gmm.weights_):
        print str(p)+'\t'+str(mu)+'\t'+str(sd)
        gmm_mean += p * mu
    print 'mix_mean\t'+str(gmm_mean)
    print 'general_mean\t'+str(mean)
    print 'AIC\t'+str(aic)
    print 'BIC\t'+str(bic)

    # prediction
    Cls = gmm.predict(np.reshape(score,(samp,1)))
    Prob = gmm.predict_proba(np.reshape(score,(samp,1)))
    Mean = gmm.means_
    Weight = gmm.weights_
    out = 'sgRNA\tGene\tFunc\tClinic\tExon\tScore\tClass\tProb\tMean\tWeight\n'
    open(outfile,'w').writelines(out)
    for j in range(0, len(score)):
        sg, ge, func, clinc, exon, sc, cls, prob = SG[j], GE[j], FUNC[j], str(CLINC[j]), str(EXON[j]), str(score[j]), str(Cls[j]), Prob[j]
        pb, mean, weight = '', '', ''
        for k in range(0, len(Mean)):
            pb += str(prob[k])+';'
            mean += str(Mean[k][0])+';'
            weight += str(Weight[k])+';'
        if clinc == 'nan':
            clinc = 'NA'
        if exon == 'nan':
            exon = 'NA' 
        out = sg+'\t'+ge+'\t'+func+'\t'+clinc+'\t'+exon+'\t'+sc+'\t'+cls+'\t'+pb+'\t'+mean+'\t'+weight+'\n'
        open(outfile,'a').writelines(out)

    # plot and compare
    fig = plt.figure(figsize = (6,3))
    fig.set_facecolor('white')
    ax = plt.subplot(111)
    #q_inf = float(pd.DataFrame(score).quantile(0.0025));
    #q_sup = float(pd.DataFrame(score).quantile(0.9975));
    q_inf = np.min(score)
    q_sup = np.max(score)
    ax.set_xlim([q_inf, q_sup])
    # empirical pdf of data
    nb = int(10*np.log(samp))
    score_f, score_n = score[score <= -0.5], score[score > -0.5]
    #ax.hist(score_f, bins = nb*2, density = False, color = '#0072B2', edgecolor = 'k', label = "Functional",  alpha=0.4)
    #ax.hist(score_n, bins = nb, density = False, color = '#999999', edgecolor = 'k', label = "Nonfunctional",  alpha=0.4)
    ax.hist(score, bins = nb, density = True, color = 'grey', edgecolor = 'k', label = "Score",  alpha=0.4)
    #sns.kdeplot(score, shade=True, color = 'grey', label = "Empirical").legend()
    # Normal fit
    #ax.plot(score_grid, normal, color = 'blue', lw = 1.0, label = "Normal fit")
    # Gaussian Mixture fit
    colors = ['orange', 'red', 'blue', 'green', 'purple', 'yellow']
    colors = ['blue', 'green', 'red', 'orange', 'yellow']
    Mean_tmp = []
    for i in range(0, num):
        Mean_tmp.append(Mean[i][0])
    Mean_tmp.sort()
    m_rank = {}
    for j in range(0, num):
        for k in range(0, num):
            if Mean_tmp[j] == Mean[k][0]:
                m_rank[j] = k
    for l in range(0, num):
        gauss_mixt_tmp = gauss_mixt[m_rank[l]]
        ax.plot(score_grid, gauss_mixt_tmp, color = colors[l], lw = 1.0, label = "GMM"+str(l+1))
    # legend
    ax.legend(loc='upper left');
    plt.tight_layout()
    #plt.show()
    plt.savefig(outfig,format='pdf',bbox_inches='tight')    


if __name__ == "__main__":
    main()
