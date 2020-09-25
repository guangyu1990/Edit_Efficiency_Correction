#!/usr/bin/python

import click
import re
import os
import math
import pandas as pd
import numpy as np
import scipy
from scipy import stats

@click.command()
@click.option("--file","-f",help="input file")

def main(file):
    with open(file) as f:
        for line in f:
            if line.startswith('sgRNA'):
                print 'sgRNA\tGene\tFunc\tClinic\tExon\tScore\tPvalue'
                continue
            inf = line.strip().split('\t')
            cls, prob, mean = inf[6], inf[7].split(';'), inf[8].split(';')
            m_value, mr_value, m_rank = [], [], {}
            num = len(mean)
            for i in range(0, num-1):
                m_value.append(float(mean[i]))
                mr_value.append(float(mean[i]))
            m_value.sort()
            for j in range(0, num-1):
                for k in range(0, len(mr_value)):
                    if m_value[j] == mr_value[k]:
                        m_rank[j] = k
            pvalue = float(prob[m_rank[0]])+float(prob[m_rank[1]])
            #pvalue = 1-(float(prob[m_rank[0]])+float(prob[m_rank[1]])+float(prob[m_rank[2]]))
            print '\t'.join(inf[0:6])+'\t'+str(pvalue)



if __name__ == "__main__":
    main()
