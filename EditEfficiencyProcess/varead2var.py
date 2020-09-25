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
@click.option("--sgfile","-s",help="sg-same file")
@click.option("--lib","-l",help="lib type")
@click.option("--outfile","-o",help="out put file")

def main(vrfile, sgfile, lib, outfile):
    global NT, SG
    NT = {'A':'G', 'T':'C', 'G':'A', 'C':'T'}
    lib = lib.split('_')[0]
    SG = read_sfile(sgfile)
    Var, Read = read_vfile(vrfile, lib)
    for vrd in Read:
        print vrd+'\t'+';'.join(Read[vrd][0:len(Read[vrd])])
    
    out = 'sgRNA\tpos\tref\talt\tdepth\tfreq\tflag\n'
    open(outfile,'w').writelines(out)
    OT = open(outfile, 'a+')
    for vs in Var:
        vs_inf = vs.split('\t')
        vs_sg, vs_pos, vs_dep, vs_ref, vs_alt, vs_var = vs_inf[0], vs_inf[1], vs_inf[2], vs_inf[3], vs_inf[4], len(Var[vs])
        
        vs_freq = float(vs_var)/float(vs_dep)
        varinf = vs_pos+'\t'+vs_ref+'\t'+vs_alt+'\t'+vs_dep+'\t'+str(vs_freq)
        tgpos, tgref, flag = '', '', ''
        if 'BRCA' in vs_sg:
            for sg_tmp in SG[vs_sg]:
                tgpos = sg_tmp.split('|')[1]
                tgref = sg_tmp.split('_')[1]
                if tgref == 'G':
                    tgref = 'C'
                elif tgref == 'T':
                    tgref = 'A'
                if tgref in lib:
                    #print vs_sg+'\t'+sg_tmp
                    vs_sg = sg_tmp
                    break                
        else:
            tgpos = vs_sg.split('|')[1]
            tgref = 'C'
        if vs_ref == tgref and NT[tgref] == vs_alt:
            if int(tgpos)+20 == int(vs_pos):
                flag = 'target'
                out = vs_sg+'\t'+varinf+'\t'+flag+'\n'
                OT.writelines(out)
            elif int(vs_pos) >= 24 and int(vs_pos) <= 28:
                flag = 'window'
                out = vs_sg+'\t'+varinf+'\t'+flag+'\n'
                OT.writelines(out)
            elif int(vs_pos) == int(tgpos):
                flag = 'seq target'
                out = vs_sg+'\t'+varinf+'\t'+flag+'\n'
                OT.writelines(out)
            elif int(vs_pos) >=4 and int(vs_pos) <= 8:
                flag = 'seq window'
                out = vs_sg+'\t'+varinf+'\t'+flag+'\n'
                OT.writelines(out)
            elif int(vs_pos) <= 20:
                flag = 'seq seed/wrong position'
                out = vs_sg+'\t'+varinf+'\t'+flag+'\n'
                OT.writelines(out)
            else:
                flag = 'seed/wrong position'
                out = vs_sg+'\t'+varinf+'\t'+flag+'\n'
                OT.writelines(out)
        else:
            flag = 'wrong base'
            out = vs_sg+'\t'+varinf+'\t'+flag+'\n'
            OT.writelines(out)
    OT.close()
    
def read_vfile(vrfile, lb):
    Vr, Rd = {}, {}
    with open(vrfile) as vf:
        for vline in vf:
            if vline.startswith('sgRNA'):
                continue
            vinf = vline.strip().split('\t')
            vsg, vpos, vdep, vref, valt, vread = vinf[0], vinf[1], vinf[2], vinf[3], vinf[4], vinf[5]
            vid = '\t'.join(vinf[0:5])
            Vr.setdefault(vid, [])
            Vr[vid].append(vread)
            if 'BRCA' in vsg:
                t_sginf = SG[vsg]
                for t_s in t_sginf:
                    tref = t_s.split('_')[1]
                    if tref == 'G':
                        tref = 'C'
                    elif tref == 'T':
                        tref = 'A'
                    if int(vpos) >= 24 and int(vpos) <= 28:
                        if vref == tref and valt == NT[tref] and tref in lb:
                            Rd.setdefault(vread, [])
                            r_var = t_s+':'+':'.join(vinf[1:5])
                            Rd[vread].append(r_var)
                    
    return Vr, Rd

def read_sfile(sfile):
    Sg = {}
    with open(sfile) as sf:
        for sline in sf:
            sinf = sline.strip().split('\t')
            s, sgp = sinf[0], sinf[1]
            Sg.setdefault(s, {})
            if ';' in sgp:
                sgp_inf = sgp.split(';')
                for sk in sgp_inf:
                    Sg[s][sk] = 1
            else:
                Sg[s][sgp] = 1
    return Sg 

if __name__ == "__main__":
    main()
