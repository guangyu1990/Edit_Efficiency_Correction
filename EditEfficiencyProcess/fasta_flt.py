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

def main(file):
    FQ = open(file)
    while 1:
        line1 = FQ.readline()
        if line1.startswith('>') and 'first' in line1:
            qinf1 = line1.strip()
            qseq1 = FQ.readline().strip()
            qinf2 = FQ.readline().strip()
            qseq2 = FQ.readline().strip()
            if len(qseq1) == 20 and len(qseq2) == 20:
                print qinf1+'\n'+qseq1+'\n'+qinf2+'\n'+qseq2
        if not line1:
            break
    FQ.close()








if __name__ == "__main__":
    main()
