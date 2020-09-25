#!/usr/bin/python

import click
import re
import os
import math
import pandas as pd
import numpy as np

@click.command()
@click.option("--file","-f",help="input file")

def main(file):
    arraylib, arrayd0r1, arrayd0r2, arrayd5r1, arrayd5r2, arrayd10r1, arrayd10r2 = [], [], [], [], [], [], []
    with open(file) as cf:
        for cline in cf:
            cinf = cline.strip().split('\t')
            if cline.startswith('sgRNA'):
                continue
            if 'Negative' not in cinf[2]:
                continue
            #print cinf[2]
            mut, lib, d0r1, d0r2, d5r1, d5r2, d10r1, d10r2 = cinf[0], cinf[5], cinf[6], cinf[7], cinf[8], cinf[9], cinf[10], cinf[11]
            arraylib.append(float(lib))
            arrayd0r1.append(float(d0r1))
            arrayd0r2.append(float(d0r2))
            arrayd5r1.append(float(d5r1))
            arrayd5r2.append(float(d5r2))
            arrayd10r1.append(float(d10r1))
            arrayd10r2.append(float(d10r2))

    sjlib, sjd0r1, sjd0r2, sjd5r1, sjd5r2, sjd10r1, sjd10r2 = 1, correct_median(arrayd0r1, arrayd0r2, arrayd5r1, arrayd5r2, arrayd10r1, arrayd10r2), correct_median(arrayd0r2, arrayd0r1, arrayd5r1, arrayd5r2, arrayd10r1, arrayd10r2), correct_median(arrayd5r1, arrayd5r2, arrayd0r1, arrayd0r2, arrayd10r1, arrayd10r2), correct_median(arrayd5r2, arrayd5r1, arrayd0r1, arrayd0r2, arrayd10r1, arrayd10r2), correct_median(arrayd10r1, arrayd10r2, arrayd0r1, arrayd0r2, arrayd5r1, arrayd5r2), correct_median(arrayd10r2, arrayd10r1, arrayd0r1, arrayd0r2, arrayd5r1, arrayd5r2)
    #print str(sjd0r1)+'\t'+str(sjd0r2)+'\t'+str(sjd10r1)+'\t'+str(sjd10r2)
    with open(file) as f:
        for line in f:
            inf = line.strip().split('\t')
            if line.startswith('sgRNA'):
                print '\t'.join(inf[0:12])+'\tlib_ratio\tday0_rep1_ratio\tday0_rep2_ratio\tday5_rep1_ratio\tday5_rep2_ratio\tday10_rep1_ratio\tday10_rep2_ratio\tday5_day0_rep1\tday5_day0_rep2\tday10_day0_rep1\tday10_day0_rep2\tday10_day5_rep1\tday10_day5_rep2\tday5_day0\tday10_day0\tday10_day5\tday5_lib\tday10_lib'
                continue
            ctlib, ctd0r1, ctd0r2, ctd5r1, ctd5r2, ctd10r1, ctd10r2 = float(inf[5]), float(inf[6]), float(inf[7]), float(inf[8]), float(inf[9]), float(inf[10]), float(inf[11])
            if ctd0r1*ctd0r2*ctd5r1*ctd5r2*ctd10r1*ctd10r2 == 0:
                continue
            rtlib, rtd0r1, rtd0r2, rtd5r1, rtd5r2, rtd10r1, rtd10r2 = float(ctlib/sjlib), float(ctd0r1/sjd0r1), float(ctd0r2/sjd0r2), float(ctd5r1/sjd5r1), float(ctd5r2/sjd5r2), float(ctd10r1/sjd10r1), float(ctd10r2/sjd10r2)
            #if rtlib*rtd0r1*rtd0r2*rtd5r1*rtd5r2*rtd10r1*rtd10r2 == 0:
            #    print str(rtlib)+'\t'+str(rtd0r1)+'\t'+str(rtd0r2)+'\t'+str(rtd5r1)+'\t'+str(rtd5r2)+'\t'+str(rtd10r1)+'\t'+str(rtd10r2)
            #    print inf[0]+'\t'+str(rtlib)+'\t'+str(rtd0r1)+'\t'+str(ctd0r2)+':'+str(sjd0r2)+'\t'+str(rtd5r1)+'\t'+str(rtd5r2)+'\t'+str(rtd10r1)+'\t'+str(rtd10r2)
            d5d0r1, d5d0r2, d10d0r1, d10d0r2, d10d5r1, d10d5r2, d5d0, d10d0, d10d5, d5lib, d10lib = float(rtd5r1)/float(rtd0r1), float(rtd5r2)/float(rtd0r2), float(rtd10r1)/float(rtd0r1), float(rtd10r2)/float(rtd0r2), float(rtd10r1)/float(rtd5r1), float(rtd10r2)/float(rtd5r2), float(rtd5r1+rtd5r2)/float(rtd0r1+rtd0r2), float(rtd10r1+rtd10r2)/float(rtd0r1+rtd0r2), float(rtd10r1+rtd10r2)/float(rtd5r1+rtd5r2), float(rtd5r1+rtd5r2)/float(2*rtlib), float(rtd10r1+rtd10r2)/float(2*rtlib)
            print '\t'.join(inf[0:12])+'\t'+str(rtlib)+'\t'+str(rtd0r1)+'\t'+str(rtd0r2)+'\t'+str(rtd5r1)+'\t'+str(rtd5r2)+'\t'+str(rtd10r1)+'\t'+str(rtd10r2)+'\t'+str(d5d0r1)+'\t'+str(d5d0r2)+'\t'+str(d10d0r1)+'\t'+str(d10d0r2)+'\t'+str(d10d5r1)+'\t'+str(d10d5r2)+'\t'+str(d5d0)+'\t'+str(d10d0)+'\t'+str(d10d5)+'\t'+str(d5lib)+'\t'+str(d10lib)

def correct_median(ary1, ary2, ary3, ary4, ary5, ary6):
    sj_array = []
    for i in range(0, len(ary1)):
        sji = ary1[i]/((ary1[i]*ary2[i]*ary3[i]*ary4[i]*ary5[i]*ary6[i])**(float(1)/6))
        sj_array.append(sji)
    sj = np.median(sj_array)
    return sj


if __name__ == "__main__":
    main()
