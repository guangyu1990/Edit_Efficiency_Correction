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
@click.option("--gpfile","-g",help="group file")
@click.option("--sgfile","-s",help="sgRNA file")
@click.option("--seq1file","-q1",help="seq of first sgRNA sequence file in fasta format")
@click.option("--seq2file","-q2",help="seq of second sgRNA sequence file in fasta format")
@click.option("--lib","-l",help="lib type")
@click.option("--outfile","-o",help="out put file")

def main(gpfile, sgfile, seq1file, seq2file, lib, outfile):
    global Sginf, PAM, GF, SEQ, GP, Var, Dep
    global treads, sgreads
    treads, sgreads = 0, 0
    Sginf, PAM = read_sgfile(sgfile)
    GF = read_gpfile(gpfile)
    SEQ, lib = {}, lib.split('_')[1]
    treads = read_seqfile(seq1file)
    read_seqfile(seq2file)
    GP = pick_group(GF, lib)
    Var, Dep = {}, {}
    out = 'sgRNA\tpos\tdepth\tref\talt\treadid\n'
    open(outfile,'w').writelines(out)
    Var_get()
    OT = open(outfile, 'a+')
    for v in Var:
        for r in Var[v]:
            out = v+'\t'+r+'\n'
            OT.writelines(out)
    OT.close()

    sg_ratio = float(sgreads)/float(treads)
    print 'total reads\tsg match reads\tsg ratio'
    print str(treads)+'\t'+str(sgreads)+'\t'+str(sg_ratio)
    for sd in Sginf:
        sgdep = 0
        if sd in Dep:
            sgdep = Dep[sd]
        print sd+'\t'+str(sgdep)


def Var_get():
    for sg in GP:
        sginf = sg.split('\t')
        if len(sginf) == 1:
            ref_seq = Sginf[sginf[0]]
            depth = len(GP[sg])
            Dep[sginf[0]] = depth
            for mrd in GP[sg]:
                aln_seq1 = SEQ[mrd]['first']
                aln_seq2 = SEQ[mrd]['second']
                #print ref_seq+'\t'+aln_seq1+'\t'+aln_seq2
                var_call(ref_seq, aln_seq1, aln_seq2, mrd, depth, sginf[0])
        else:
            GP_TMP = group_split(sginf, GP[sg])             

def group_split(sgf, gprd):
    Var_tmp = {}
    global sgreads
    for gr in gprd:
        al_sq1_tmp = SEQ[gr]['first']
        al_sq2_tmp = SEQ[gr]['second']
        al_pam_tmp = gr.split('_')[1]
        sg_vnum = []
        V_tmp = {}
        for s in sgf:
            rf_sq_tmp, rf_pam_tmp = Sginf[s], PAM[s]
            v_num = 0
            if al_pam_tmp != rf_pam_tmp:
                v_num = 10000
            for j in range(0, len(rf_sq_tmp)):
                vid_tmp = ''
                if al_sq1_tmp[j] != rf_sq_tmp[j]:
                    v_num += 1
                    vid_tmp = s+'\t'+str(j+1)+'\t'+rf_sq_tmp[j]+'\t'+al_sq1_tmp[j]
                    V_tmp.setdefault(s, [])
                    V_tmp[s].append(vid_tmp)
                if al_sq2_tmp[j] != rf_sq_tmp[j]:
                    vid_tmp = s+'\t'+str(j+1+20)+'\t'+rf_sq_tmp[j]+'\t'+al_sq2_tmp[j]
                    V_tmp.setdefault(s, [])
                    V_tmp[s].append(vid_tmp)
            sg_vnum.append(v_num)
        if min(sg_vnum) < 10000:
            sgreads += 1
            idx_ok = sg_vnum.index(min(sg_vnum))
            Dep.setdefault(sgf[idx_ok], 0)
            Dep[sgf[idx_ok]] += 1
            Var_tmp.setdefault(sgf[idx_ok], {})
            if sgf[idx_ok] in V_tmp:
                for vtp in V_tmp[sgf[idx_ok]]:
                    Var_tmp[sgf[idx_ok]].setdefault(vtp, [])
                    Var_tmp[sgf[idx_ok]][vtp].append(gr)
    for l in Var_tmp:
        for k in Var_tmp[l]:
            k_inf = k.split('\t')
            vid_ok = '\t'.join(k_inf[0:2])+'\t'+str(Dep[l])+'\t'+'\t'.join(k_inf[2:len(k_inf)])
            Var.setdefault(vid_ok, [])
            Var[vid_ok] = Var_tmp[l][k]

def var_call(rf_sq, al_sq1, al_sq2, mr, dep, sgid):
    for i in range(0, len(rf_sq)):
        vid = ''
        if al_sq1[i] != rf_sq[i]:
            vid = sgid+'\t'+str(i+1)+'\t'+str(dep)+'\t'+rf_sq[i]+'\t'+al_sq1[i]
            Var.setdefault(vid, [])
            Var[vid].append(mr)
        if al_sq2[i] != rf_sq[i]:
            vid = sgid+'\t'+str(i+1+20)+'\t'+str(dep)+'\t'+rf_sq[i]+'\t'+al_sq2[i]
            Var.setdefault(vid, [])
            Var[vid].append(mr)

def pick_group(Gf, lb):
    Gp = {}
    global sgreads
    for cid in Gf:
        gflg = 0
        sgroup = []
        reads = {}
        for sid in Gf[cid]:
            if sid in Sginf:
                if lb == 'NG' and '_NG_' in sid:
                    gflg += 1
                    sgroup.append(sid)
                elif lb == 'NGG' and '_NG_' not in sid:
                    gflg += 1
                    sgroup.append(sid)
                elif lb == 'lib':
                    gflg += 1
                    sgroup.append(sid)
            elif '|' not in sid:
                rinf = sid.split('_')
                rid, l_b, tg = rinf[0], rinf[1], rinf[2]
                reads[rid+'_'+l_b] = 1
        if gflg > 0:
            sgref = '\t'.join(sgroup[0:len(sgroup)])
            sg_pam = ''
            if gflg == 1:
                sg_pam = PAM[sgref]
            if reads:
                for rd in reads:
                    rd_pam = rd.split('_')[1]
                    if gflg == 1:
                        if sg_pam == rd_pam:
                            sgreads += 1
                            Gp.setdefault(sgref, [])
                            Gp[sgref].append(rd)
                    else:
                        Gp.setdefault(sgref, [])
                        Gp[sgref].append(rd)
    return Gp

def read_sgfile(sfile):
    Sg, Pam = {}, {}
    with open(sfile) as sf:
        for sline in sf:
            sinf = sline.strip().split('\t')
            sg, spos, seq = sinf[0], sinf[1], sinf[2]
            tag, pam = seq[21:41], seq[146:149]
            Sg[sg+'|'+spos] = tag
            Pam[sg+'|'+spos] = pam
    return Sg, Pam

def read_gpfile(gfile):
    Gf = {}
    cltid = ''
    with open(gfile) as gf:
        for gline in gf:
            if gline.startswith('>'):
                cltid = gline.strip()
                Gf.setdefault(cltid, [])
            else:
                ginf = gline.strip().split('>')[1].split('...')[0]
                Gf[cltid].append(ginf)
    return Gf

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
            qid, qlb, qtg = qinf[0], qinf[1], qinf[2]
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
