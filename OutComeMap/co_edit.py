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
@click.option("--vrfile","-v",help="var file")
@click.option("--sgfile","-s",help="sgRNA seq file")
@click.option("--msfile","-m",help="match sgRNA file")
@click.option("--outfile","-o",help="out put file")

def main(vrfile, sgfile, msfile, outfile):
    global NT
    NT = {'A':'A', 'T':'A', 'G':'C', 'C':'C'}
    SG = read_sfile(sgfile)
    MS = read_mfile(msfile)
    Var = read_vfile(vrfile)
    out = 'sgRNA\tseq\tdepth\tvar_pos\ttarget\twindow\tco_edit\ttotal\n'
    open(outfile,'w').writelines(out)
    OT = open(outfile, 'a+')
    for v in Var:
        total_var = {}
        dep = Var[v]['depth']
        del Var[v]['depth']
        seq = ''
        for m in MS[v]:
            if m in SG:
                seq = SG[m]
        var_window = get_window(v, seq)
        t_pos, t_alt = v.split('|')[1], ''
        w_alt, co_edit = '', ''
        for k in Var[v]:
            total_var[k+'|'] = Var[v][k]
            if 'co' in k:
                co_edit += k+'|'+str(Var[v][k])+';'
            else:
                if k == t_pos:
                    t_alt = str(Var[v][k])
                else:
                    w_alt += k+'|'+str(Var[v][k])+';' 
        if t_alt == '':
            t_alt = 'NA'
        if w_alt == '':
            w_alt = 'NA'
        total_out = ''
        for t in sorted(total_var.items(), key=lambda x: x[1], reverse=True):
            total_out += t[0]+str(t[1])+';'
        out = v+'\t'+seq+'\t'+dep+'\t'+':'.join(var_window)+'\t'+t_alt+'\t'+w_alt+'\t'+co_edit+'\t'+total_out+'\n'   
        OT.writelines(out)
    OT.close()     

def get_window(id, sq):
    wd = []
    s = NT[id.split('_')[1]]
    for i in range(3, 8):
        if sq[i] == s:
            wd.append(str(i+1))
    return wd

def read_vfile(vrfile):
    Vr = {}
    with open(vrfile) as vf:
        for vline in vf:
            vinf = vline.strip().split('\t')
            vread, vid, vdp, vtps, vsg = vinf[0], vinf[1].split(':')[0], vinf[1].split(':')[2], vinf[1].split('|')[1], []
            Vr.setdefault(vid, {})
            Vr[vid]['depth'] = vdp
            if ';' in vinf[1]:
                vsg = vinf[1].split(';')
            else:
                vsg.append(vinf[1])
            co_pos = []
            for var in vsg:
                vpos = int(var.split(':')[1])-20
                Vr[vid].setdefault(str(vpos), 0)
                Vr[vid][str(vpos)] += 1
                co_pos.append(vpos)
            co_pos.sort()
            co_pos = map(str, co_pos)
            co_pos = 'co:'+':'.join(co_pos)
            #print co_pos
            Vr[vid].setdefault(co_pos, 0)
            Vr[vid][co_pos] += 1

    return Vr

def read_sfile(sfile):
    Sg = {}
    FQ = open(sfile)
    while 1:
        sline = FQ.readline()
        if sline.startswith('>'):
            sgid = sline.strip().replace('>','')
            sgseq = FQ.readline().strip()
            Sg[sgid] = sgseq
        if not sline:
            break
    FQ.close()
    return Sg 

def read_mfile(mfile):
    Ms = {}
    with open(mfile) as mf:
        for mline in mf:
            minf = mline.strip().split('\t')
            m, mgp = minf[0], minf[1]
            Ms.setdefault(m, {})
            if ';' in mgp:
                mgp_inf = mgp.split(';')
                for mk in mgp_inf:
                    Ms[m][mk] = 1
            else:
                Ms[m][mgp] = 1
    return Ms


if __name__ == "__main__":
    main()
