#!/usr/bin/python

import click
import re
import os
import math
import pandas as pd
import numpy as np
import scipy as sp
from multiprocessing import Process,Manager

@click.command()
@click.option("--gpfile","-g",help="group file")
@click.option("--sgfile","-s",help="sgRNA file")
@click.option("--seqfile","-q",help="seq of sgRNA sequence file in fasta format")
@click.option("--outfile","-o",help="out put file")

def main(gpfile, sgfile, seqfile, outfile):
    global Sginf, GF, SEQ, GP, Var, Dep
    Sginf = read_sgfile(sgfile)
    GF, Ngroup = read_gpfile(gpfile)
    SEQ, Lib = read_seqfile(seqfile)
    GP, sggp, usggp, treads, mreads = pick_group(GF, Lib)
    print 'total group\tsgRNA group\tuniq sgRNA group\ttotal reads\tmatch reads'
    print str(Ngroup)+'\t'+str(sggp)+'\t'+str(usggp)+'\t'+str(treads)+'\t'+str(mreads)
    Var, Dep = {}, {}
    out = 'sgRNA\tpos\tdepth\tref\talt\treadid\n'
    open(outfile,'w').writelines(out)
    Var_get()
    #for v in Var:
    #    for r in Var[v]:
    #        out = v+'\t'+r+'\n'
    #        open(outfile,'a').writelines(out)
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
            #for mrd in GP[sg]:
            #    aln_seq1 = SEQ[mrd]['first']
            #    aln_seq2 = SEQ[mrd]['second']
            #    #print ref_seq+'\t'+aln_seq1+'\t'+aln_seq2
            #    var_call(ref_seq, aln_seq1, aln_seq2, mrd, depth, sginf[0])
        else:
            GP_TMP = group_split(sginf, GP[sg])             

def group_split(sgf, gprd):
    Var_tmp = {}
    for gr in gprd:
        al_sq1_tmp = SEQ[gr]['first']
        al_sq2_tmp = SEQ[gr]['second']
        sg_vnum = []
        V_tmp = {}
        for s in sgf:
            rf_sq_tmp = Sginf[s]
            v_num = 0
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
        idx_ok = sg_vnum.index(min(sg_vnum))
        Dep.setdefault(sgf[idx_ok], 0)
        Dep[sgf[idx_ok]] += 1
        #Var_tmp.setdefault(sgf[idx_ok], {})
        #if sgf[idx_ok] in V_tmp:
        #    for vtp in V_tmp[sgf[idx_ok]]:
        #        Var_tmp[sgf[idx_ok]].setdefault(vtp, [])
        #        Var_tmp[sgf[idx_ok]][vtp].append(gr)
    #for l in Var_tmp:
    #    for k in Var_tmp[l]:
    #        k_inf = k.split('\t')
    #        vid_ok = '\t'.join(k_inf[0:2])+'\t'+str(Dep[l])+'\t'+'\t'.join(k_inf[2:len(k_inf)])
    #        Var.setdefault(vid_ok, [])
    #        Var[vid_ok] = Var_tmp[l][k]

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

def pick_group(Gf, lib):
    sgp, ugp, tnum, mnum = 0, 0, 0, 0
    Gp = {}
    for cid in Gf:
        gflg = 0
        sgroup = []
        reads = {}
        for sid in Gf[cid]:
            if sid in Sginf:
                if lib == 'NG' and '_NG_' in sid:
                    gflg += 1
                    sgroup.append(sid)
                elif lib == 'NGG' and '_NG_' not in sid:
                    gflg += 1
                    sgroup.append(sid)
            elif '|' not in sid:
                tnum += 1
                rinf = sid.split('_')
                rid, tg, lb = rinf[0], rinf[1], rinf[2]
                if rid+'_'+lb in reads:
                    reads[rid+'_'+lb] += 1
                else:
                    reads[rid+'_'+lb] = 1
        if gflg > 0:
            sgp += 1
            sgref = '\t'.join(sgroup[0:len(sgroup)])
            if reads:
                for rd in reads:
                    if reads[rd] > 1:
                        mnum += 1
                        Gp.setdefault(sgref, [])
                        Gp[sgref].append(rd)
        if gflg == 1:
            ugp += 1
    return Gp, sgp, ugp, tnum, mnum

def read_sgfile(sfile):
    Sg = {}
    with open(sfile) as sf:
        for sline in sf:
            sinf = sline.strip().split('\t')
            sg, spos, seq = sinf[0], sinf[1], sinf[2]
            tag = seq[21:41]
            Sg[sg+'|'+spos] = tag
    return Sg

def read_gpfile(gfile):
    Gf = {}
    cltid = ''
    gp_num = 0
    with open(gfile) as gf:
        for gline in gf:
            if gline.startswith('>'):
                gp_num += 1
                cltid = gline.strip()
                Gf.setdefault(cltid, [])
            else:
                ginf = gline.strip().split('>')[1].split('...')[0]
                Gf[cltid].append(ginf)
    return Gf, gp_num           

def read_seqfile(qfile):
    Sq, splb = {}, ''
    FQ = open(qfile)
    while 1:
        line1 = FQ.readline()
        if line1.startswith('>'):
            qinf = line1.strip().split('_')
            qseq = FQ.readline().strip()
            qid, qtg, qlb = qinf[0], qinf[1], qinf[2]
            splb = qlb
            sqid = qid.replace('>','')+'_'+qlb
            Sq.setdefault(sqid, {})
            Sq[sqid][qtg] = qseq
        if not line1:
            break
    FQ.close()
    return Sq, splb

if __name__ == "__main__":
    main()
