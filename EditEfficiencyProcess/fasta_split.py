#!/usr/bin/python

import click
import re
import os
import math
import pandas as pd
import numpy as np
import scipy as sp

@click.command()
@click.option("--file","-f",help="input fasta file")
@click.option("--out","-o",help="output fasta file")

def main(file, out):
    outfile1 = out+'.first.fasta'
    outfile2 = out+'.second.fasta'
    out1 = ''
    out2 = ''
    open(outfile1,'w').writelines(out1)
    open(outfile2,'w').writelines(out2)
    FQ = open(file)
    O1 = open(outfile1, 'a+')
    O2 = open(outfile2, 'a+')
    while 1:
        line1 = FQ.readline()
        if line1.startswith('>') and 'first' in line1:
            qinf1 = line1.strip()
            qseq1 = FQ.readline().strip()
            qinf2 = FQ.readline().strip()
            qseq2 = FQ.readline().strip()
            out1 = qinf1+'\n'+qseq1+'\n'
            out2 = qinf2+'\n'+qseq2+'\n'
            #open(outfile1,'a').writelines(out1)
            #open(outfile2,'a').writelines(out2)
            O1.writelines(out1)
            O2.writelines(out2)
        if not line1:
            break
    FQ.close()
    O1.close()
    O2.close()



if __name__ == "__main__":
    main()
