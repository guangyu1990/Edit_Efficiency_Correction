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

def main(vrfile):
    with open(vrfile) as vf:
        for vline in vf:
            if vline.startswith('sgRNA'):
                print 'sgRNA\tseq\tdepth\tvar_pos\tpick_pos\tpick_aa\tpick_type\tpick_support\ttarget\twindow\tco_edit\ttotal'
                continue
            vinf = vline.strip().split('\t')
            vsg, vvar = vinf[0], vinf[7]
            sg, tpos = vsg.split('_sg1')[0], vsg.split('|')[1]
            var_inf = vvar.split(';')
            AA, AF = {}, {}
            if 'UTR' in vvar:
                uaa = re.search(r'UTR\d', vvar).group()
                PS, PS_rd = {}, {}
                for u in range(0, len(var_inf)-1):
                    utr_inf = var_inf[u].split('|')
                    upos, ualt = utr_inf[0].replace('co:',''), utr_inf[1]
                    for p in upos.split(':'):
                        PS.setdefault(p, 0)
                        PS[p] += float(ualt)
                        PS_rd.setdefault(p, [])
                        PS_rd[p].append('co:'+upos)
                u_maxp, u_max = sorted(PS.items(), key = lambda kv:(kv[1], kv[0]), reverse=True)[0][0], sorted(PS.items(), key = lambda kv:(kv[1], kv[0]), reverse=True)[0][1]
                umaxp = '|'.join(PS_rd[u_maxp])
                out = '\t'.join(vinf[0:4])+'\t'+umaxp+'\t'+uaa+'\t'+uaa+'\t'+str(u_max)+'\t'+'\t'.join(vinf[4:len(vinf)])
                print out
                continue
            for i in range(0, len(var_inf)-1):
                all_inf = var_inf[i].split('|')
                pos, alt, nt, aa, mtyp, aa_fix = all_inf[0], all_inf[1], all_inf[2], all_inf[3].split(':'), all_inf[4].split(':'), []
                for a_tmp in aa:
                    a_fix = ''
                    if not re.search(r'([A-Z]\d+)', a_tmp):
                        a_fix = a_tmp
                    else:
                        a_fix = re.search(r'([A-Z]\d+)', a_tmp).group()
                    aa_fix.append(a_fix)
                    AA.setdefault(a_fix, {})
                    AA[a_fix].setdefault(a_tmp, {})
                    AA[a_fix][a_tmp].setdefault('num', 0)
                    AA[a_fix][a_tmp]['num'] += float(alt)
                    AA[a_fix][a_tmp].setdefault('pos', [])
                    AA[a_fix][a_tmp]['pos'].append(pos)
                    AA[a_fix][a_tmp]['mtype'] = get_aatype(a_tmp)
                aa_fix = list(set(aa_fix))
                for j in range(0, len(aa_fix)):
                    AF.setdefault(aa_fix[j], 0)
                    AF[aa_fix[j]] += float(alt)
                
            alt_max, aa_max, pos_max, mtyp_max = 0, '', '', ''
            for a in AF:
                if AF[a] > alt_max:
                    alt_max = AF[a]
                    aa_max_ary, aa_maxnum_ary, pos_max_ary, mtyp_max_ary = [], {}, [], {}
                    for al in AA[a]:
                        aa_maxnum_ary[al] = AA[a][al]['num']
                        pos_max_ary += AA[a][al]['pos']
                        mtyp_max_ary.setdefault(AA[a][al]['mtype'], 0)
                        mtyp_max_ary[AA[a][al]['mtype']] += AA[a][al]['num']
                    for aal in sorted(aa_maxnum_ary.items(), key = lambda kv:(kv[1], kv[0]), reverse=True):
                        aa_max_ary.append(aal[0])
                    #print mtyp_max_ary
                    pos_max_ary = list(set(pos_max_ary))
                    aa_max = '|'.join(aa_max_ary)
                    pos_max = '|'.join(pos_max_ary)
                    mtyp_max = sorted(mtyp_max_ary.items(), key = lambda kv:(kv[1], kv[0]), reverse=True)[0][0]
            #print mtyp_max
            out = '\t'.join(vinf[0:4])+'\t'+pos_max+'\t'+aa_max+'\t'+mtyp_max+'\t'+str(alt_max)+'\t'+'\t'.join(vinf[4:len(vinf)])
            print out

def get_aatype(aa_change):
    mut_type = 'NA'
    if re.search(r'([A-Z]\d+)', aa_change):
        raw = re.findall(r'[A-Z]', aa_change)[0]
        new = re.findall(r'[A-Z]', aa_change)[1]
        if new == raw:
            mut_type = 'synonymous SNV'
        elif new == 'X':
            mut_type = 'stopgain'
        elif raw == 'X':
            mut_type = 'stoploss'
        else:
            mut_type = 'nonsynonymous SNV'
    else:
        mut_type = aa_change
    return mut_type

def merge_mytype(marry):
    maa = 'NA'
    if 'stopgain' in marry:
        maa = 'stopgain'
    elif 'stoploss' in marry:
        maa = 'stoploss'
    elif 'splicing' in marry:
        maa = 'splicing'
    elif 'nonsynonymous SNV' in marry:
        maa = 'nonsynonymous SNV'
    elif 'UTR5' in marry:
        maa = 'UTR5'
    elif 'UTR3' in marry:
        maa = 'UTR3'
    elif 'synonymous SNV' in marry:
        maa = 'synonymous SNV'
    else:
        maa = marry[0]
    return maa

if __name__ == "__main__":
    main()
