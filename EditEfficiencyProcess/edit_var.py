#!/usr/bin/python

import os
import gzip
import click

@click.command()
@click.option("--fastq","-fq",help="merge fastq")
@click.option("--sgfile","-s",help="sgRNA sequence file")
@click.option("--rbfile","-r",help="recombine file from matched lib control")
@click.option("--outfile","-o",help="out file")

def main(sgfile, fastq, rbfile, outfile):
    Sginf, anchor = read_sgfile(sgfile)
    #print anchor
    Dseqs, treads = read_fastq(fastq)
    Rseqs, rreads = read_rbfile(rbfile)
    Stat, Var = {}, {}
    Recobine = {}
    ach_reads, recb_reads = 0, 0

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
                sg_merg = sg_tmp1 + sg_tmp2
                if sg_merg in Rseqs:
                    recb_reads += 1
                    break
                else:
                    v_rge = len(sg_tmp1)
                    if pam_tmp == 'GG':
                        sg_tmp1 = sg_tmp1+'_NGG'
                    if sg_tmp1 in Sginf:
                        #print 'right\t'+tg
                        #print 'rigth\t'+sequence
                        idx = Sginf[sg_tmp1]
                        #print idx+'\t'+sg_tmp1+'\t'+sg_tmp2
                        if idx in Stat:
                            Stat[idx] += 1
                        else:
                            Stat[idx] = 1
                        for i in range(0, v_rge):
                            if sg_tmp2[i] != sg_tmp1[i]:
                                var_idx = Sginf[sg_tmp1]+'\t'+str(i+1)+'\t'+sg_tmp1[i]+'\t'+sg_tmp2[i]
                                if var_idx in Var:
                                    Var[var_idx] += 1
                                else:
                                    Var[var_idx] = 1
                    #else:
                        #print 'wrong\t'+tg
                        #print 'wrong\t'+sequence
                break
    tsg_num = 0
    for d in Sginf:
        tmp = Sginf[d]
        tot = 0
        if tmp in Stat:
            tot = Stat[tmp]
        tsg_num += tot
        out = tmp+' total_number:\t'+str(tot)
        print out
    achor_rate = float(ach_reads)/float(treads)
    tsg_rate = float(tsg_num)/float(ach_reads)
    recb_rate = float(recb_reads)/float(ach_reads)
    print '##total reads\tachor reads\tachor rate\ttotal sgRNA reads\ttotal sgRNA rate\trecb reads\trecb rate'
    print str(treads)+'\t'+str(ach_reads)+'\t'+str(achor_rate)+'\t'+str(tsg_num)+'\t'+str(tsg_rate)+'\t'+str(recb_reads)+'\t'+str(recb_rate)

    var_out = '#sgRNA\tpos\tref\talt\tdepth\tvars\tfreq\n'
    open(outfile,'w').writelines(var_out)
    for v in Var:
        vinf = v.split('\t')
        vsg, vpos, vref, valt = vinf[0], vinf[1], vinf[2], vinf[3]
        vdep, vvars = Stat[vsg], Var[v]
        vfreq = float(vvars)/float(vdep)
        var_out = vsg+'\t'+vpos+'\t'+vref+'\t'+valt+'\t'+str(vdep)+'\t'+str(vvars)+'\t'+str(vfreq)+'\n'
        open(outfile,'a').writelines(var_out)

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
    Dqs = {}
    FQ = gzip.open(qfile)
    while 1:
        line1 = FQ.readline()
        if line1.startswith('@'):
            flag = line1.strip()
            seq = FQ.readline().strip()
            Dqs[flag] = seq
            trds += 1
        if not line1:
            break
    FQ.close()
    return Dqs, trds

def read_rbfile(rfile):
    Rqs = {}
    rrds = 0
    with open(rfile) as rf:
        for rline in rf:
            if rline.startswith('#'):
                continue
            rinf = rline.strip().split('\t')
            rseq, rrate = rinf[0], rinf[2]
            Rqs[rseq] = rrate
            rrds += 1
    return Rqs, rrds

if __name__ == "__main__":
    main()
