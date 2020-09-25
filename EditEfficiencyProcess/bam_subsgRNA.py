#!/usr/bin/python

import click
import re
import os
import math
import pandas as pd
import numpy as np
import scipy as sp
import pysam
from multiprocessing import Process,Manager

@click.command()
@click.option("--bam","-b",help="bam file")
@click.option("--out","-o",help="out put file")

def main(bam, out):
    bamfile = pysam.AlignmentFile(bam, "rb")
    treads, mreads, mlreads = 0, 0, 0
    out1, out2 = '', ''
    outfile1, outfile2 = out+'_sgseq-first.fasta', out+'_sgseq-second.fasta'
    open(outfile1,'w').writelines(out1)
    open(outfile2,'w').writelines(out2)
    OT1 = open(outfile1, 'a+')
    OT2 = open(outfile2, 'a+')
    for read in bamfile:
        treads += 1
        name = read.query_name
        flag = int(read.flag)
        ref = int(read.reference_id)
        pos = int(read.reference_start)
        mq = int(read.mapping_quality)
        cigar = read.cigarstring
        seq = read.query_sequence
        #qua = np.array(read.query_qualities)+33
        #qua_new = ''
        #for q in qua:
        #    qua_new += chr(q)
        #print name+'\t'+str(ref)+'\t'+str(pos)+'\t'+str(mq)
        if ref == 0 and pos == 0 and mq == 60 and flag == 0:
            tag_nm = int(read.get_tag('NM'))
            #print name+'\t'+str(ref)+'\t'+str(pos)+'\t'+str(mq)+'\t'+str(tag_nm)
            if tag_nm > 3:
                continue
            mreads += 1
            st1, ed1, st2, ed2 = get_pos(cigar, len(seq))
            pam = seq[ed2:ed2+3]
            #if pam == 'GG':
            #    pam = 'NGG'
            #else:
            #    pam = 'NG'
            #if lib != 'lib' and pam != lib:
            #    #print pam+'\t'+lib
            #    continue
            #mlreads += 1
            #print '@'+name+'\n'+seq+'\n+\n'+qua_new
            sg1 = seq[st1:ed1]
            sg2 = seq[st2:ed2]
            out1 = '>'+name+'_first_'+pam+'\n'+sg1+'\n'
            out2 = '>'+name+'_second_'+pam+'\n'+sg2+'\n'
            OT1.writelines(out1)
            OT2.writelines(out2)
    OT1.close()
    OT2.close()
    mrate = round(float(mreads)/float(treads), 4)
    #mlrate = round(float(mlreads)/float(treads), 4)
    print '##total reads\tachor match reads\tachor match rate'
    print str(treads)+'\t'+str(mreads)+'\t'+str(mrate)

def get_pos(cg, length):
    s1, e1, s2, e2 = 0, 0, 0, 0
    ach = 82
    delnum, insnum = re.findall(r'(\d+)D',cg), re.findall(r'(\d+)I',cg)
    s = re.search(r'(\d+)[A-Z]+',cg).group(1)
    for d in delnum:
        ach = ach - int(d)
    for i in insnum:
        ach = ach + int(i)
    s1, e1 = int(s)-20, int(s)
    s2, e2 = int(s)+ach+3, int(s)+ach+23
    return s1, e1, s2, e2








if __name__ == "__main__":
    main()
