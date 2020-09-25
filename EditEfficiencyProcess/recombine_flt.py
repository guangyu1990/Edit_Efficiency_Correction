#!/usr/bin/python

import click
import re
import os
import math
import pandas as pd
import numpy as np
import scipy as sp
import Levenshtein
from multiprocessing import Process,Manager

@click.command()
@click.option("--seq1file","-q1",help="seq of first sgRNA sequence file in fasta format")
@click.option("--seq2file","-q2",help="seq of second sgRNA sequence file in fasta format")
@click.option("--ident","-i",type=float,help="identity")
@click.option("--dist","-d",help="distance method",default='ratio')
@click.option("--outpre","-o",help="out put file")

def main(seq1file, seq2file, ident, dist, outpre):
    global SEQ, treads, mreads
    treads, mreads = 0, 0
    SEQ = {}
    treads = read_seqfile(seq1file)
    read_seqfile(seq2file)
    outfile1, outfile2 = outpre+'_sgseq-first_fltrb.fasta', outpre+'_sgseq-second_fltrb.fasta'
    out = ''
    open(outfile1,'w').writelines(out)
    open(outfile2,'w').writelines(out)
    OT1 = open(outfile1, 'a+')
    OT2 = open(outfile2, 'a+')
    for rd in SEQ:
        rd_first = SEQ[rd]['first']
        rd_second = SEQ[rd]['second']
        if len(rd_first) != 20 or len(rd_second) != 20:
            continue
        rd_mc_ratio = 0
        if dist == 'jaro':
            rd_mc_ratio = Levenshtein.jaro(rd_first, rd_second)
        else:
            rd_mc_ratio = Levenshtein.ratio(rd_first, rd_second)
        #print rd_mc_ratio
        if rd_mc_ratio >= ident:
            mreads += 1
            out1 = '>'+rd+'_first\n'+rd_first+'\n'
            out2 = '>'+rd+'_second\n'+rd_second+'\n'
            OT1.writelines(out1)
            OT2.writelines(out2)
    OT1.close()
    OT2.close()

    match_ratio = float(mreads)/float(treads)
    reco_ratio = float(treads-mreads)/float(treads)
    print 'total reads\tmatch reads\tmatch ratio\trecombine ratio'
    print str(treads)+'\t'+str(mreads)+'\t'+str(match_ratio)+'\t'+str(reco_ratio)

def read_seqfile(qfile):
    FQ = open(qfile)
    tnum = 0
    while 1:
        line1 = FQ.readline()
        if line1.startswith('>'):
            if 'first' in line1.strip():
                tnum += 1
            qinf = line1.strip().split('_')
            qseq = FQ.readline().strip()
            qid, qtg, qlb = qinf[0], qinf[1], qinf[2]
            sqid = qid.replace('>','')+'_'+qlb
            SEQ.setdefault(sqid, {})
            SEQ[sqid][qtg] = qseq
        if not line1:
            break
    FQ.close()
    if tnum > 0:
        return tnum

if __name__ == "__main__":
    main()
