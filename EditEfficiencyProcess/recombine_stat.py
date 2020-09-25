#!/usr/bin/python

import os
import gzip
import click

@click.command()
@click.option("--fastq","-fq",help="merge fastq")
@click.option("--sgfile","-s",help="sgRNA sequence file")
@click.option("--outfile","-o",help="out file")

def main(sgfile, fastq, outfile):
    Sginf, anchor = read_sgfile(sgfile)
    #print anchor
    Dseqs, treads = read_fastq(fastq)
    Flag, Stat = {}, {}
    Recobine = {}
    ach_reads = 0

    for tg in Dseqs:
        Tmp = {}
        sequence = Dseqs[tg]
        #print sequence
        for i in range(0, len(sequence)-82):
            seq_tmp = sequence[i:i+82]
            if seq_tmp == anchor:
                #print tg
                ach_reads += 1
                sg_tmp1 = sequence[i-20:i]
                sg_tmp2 = sequence[i+85:i+105]
                pam_tmp = sequence[i+106:i+108]
                if sg_tmp1 == sg_tmp2:
                    if pam_tmp == 'GG':
                        sg_tmp1 = sg_tmp1+'_NGG'
                    if sg_tmp1 in Sginf:
                        idx = Sginf[sg_tmp1]+'_'+sg_tmp1
                        if idx in Stat:
                            Stat[idx] += 1
                        else:
                            Stat[idx] = 1
                        if idx in Flag:
                            Flag[idx] += 1
                        else:
                            Flag[idx] = 1
                else:
                    recb_seq = sg_tmp1+sg_tmp2
                    if recb_seq in Recobine:
                        Recobine[recb_seq] += 1
                    else:
                        Recobine[recb_seq] = 1
                    if pam_tmp == 'GG':
                        sg_tmp2 = sg_tmp2+'_NGG'
                    if sg_tmp1 in Sginf:
                        idx = Sginf[sg_tmp1]+'_'+sg_tmp1
                        if idx in Stat:
                            Stat[idx] += 1
                        else:
                            Stat[idx] = 1
                    elif sg_tmp2 in Sginf:
                        idx = Sginf[sg_tmp2]+'_'+sg_tmp2
                        if idx in Stat:
                            Stat[idx] += 1
                        else:
                            Stat[idx] = 1
                break
    
    tmat = 0
    for d in Sginf:
        tmp = Sginf[d]+'_'+d
        tot, mat = 0, 0
        if tmp in Stat:
            tot = Stat[tmp]
        if tmp in Flag:
            mat = Flag[tmp]
        tmat += mat
        out = tmp+' total_number|match_number:\t'+str(tot)+'\t'+str(mat)
        print out
    achor_rate = float(ach_reads)/float(treads)
    match_rate = float(tmat)/float(ach_reads)
    print '##total reads\tachor reads\tachor rate\tmatch reads\tmatch rate'
    print str(treads)+'\t'+str(ach_reads)+'\t'+str(achor_rate)+'\t'+str(tmat)+'\t'+str(match_rate)

    rout = '##recobine seq\trecobine reads\trecobine rate\n'
    open(outfile,'w').writelines(rout)
    for rcb in Recobine:
        rcb_rate = round(float(Recobine[rcb])/float(ach_reads-tmat), 4)
        rout = rcb+'\t'+str(Recobine[rcb])+'\t'+str(rcb_rate)+'\n'
        open(outfile,'a').writelines(rout)

def read_sgfile(sfile):
    Sg, ach = {}, ''
    with open(sfile) as sf:
        for sline in sf:
            sinf = sline.strip().split('\t')
            sg, seq = sinf[0], sinf[2]
            tag = seq[21:41]
            pam = seq[147:149]
            if pam == 'GG':
                tag = tag + '_NGG'
            ach = seq[41:123].upper()
            Sg[tag] = sg
    return Sg, ach

def read_fastq(qfile):
    trds = 0
    Dseqs = {}
    FQ = gzip.open(qfile)
    while 1:
        line1 = FQ.readline()
        if line1.startswith('@'):
            flag = line1.strip()
            seq = FQ.readline().strip()
            Dseqs[flag] = seq
            trds += 1
        if not line1:
            break
    FQ.close()
    return Dseqs, trds

if __name__ == "__main__":
    main()
