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
@click.option("--anfile","-a",help="anno file")

def main(vrfile, anfile):
    Anno = read_anfile(anfile)
    ALT = {'A':'G','T':'C','G':'A','C':'T'}
    with open(vrfile) as vf:
        for vline in vf:
            if vline.startswith('sgRNA'):
                print vline.strip()
                continue
            vinf = vline.strip().split('\t')
            vsg, vvar = vinf[0], vinf[7]
            sg, tpos = vsg.split('_sg1')[0], vsg.split('|')[1]
            var_inf = vvar.split(';')
            mtyp_merge, ano_inf = Anno[sg]['mtype'], Anno[sg]['anno']
            var_out = '\t'.join(vinf[0:7])+'\t'
            for j in range(0,len(var_inf)-1):
                if 'co:' not in var_inf[j]:
                    continue
                w_pos = var_inf[j].split('|')[0]
                w_pos_tmp = w_pos.replace('co:', '')
                flg = 0
                for m in range(0, len(ano_inf)):
                    wan_pos = ano_inf[m].split('|')[0]
                    if w_pos_tmp == wan_pos:
                        flg += 1
                        var_out += var_inf[j]+'|'+'|'.join(ano_inf[m].split('|')[1:4])+';'
                if flg == 0:
                    #print sg+'\t'+w_pos
                    w_pos_co = var_inf[j].split('|')[0].replace('co:','').split(':')
                    #print w_pos_co
                    w_c, w_a, w_mt, w_p = [], [], [], []
                    for w in range(0,len(w_pos_co)):
                        for x in range(0, len(ano_inf)):
                            if w_pos_co[w] == ano_inf[x].split('|')[0]:
                                flg += 1
                                #print sg+'\t'+ano_inf[x].split('|')[0]
                                w_c.append(ano_inf[x].split('|')[1])
                                w_a.append(ano_inf[x].split('|')[2])
                                w_mt.append(ano_inf[x].split('|')[3])
                                w_p.append(w_pos_co[w])
                    if len(list(set(w_a))) == 1:
                        w_a = list(set(w_a))
                        w_mt = list(set(w_mt))
                    AA = {}
                    for y in range(0, len(w_a)):
                        aa = w_a[y]
                        if 'p.' in aa:
                            aa = re.search(r'([A-Z]\d+)',aa).group()
                        AA.setdefault(aa,[])
                        AA[aa].append(w_c[y]+'\t'+w_a[y]+'\t'+w_p[y]+'\t'+w_mt[y])
                    w_a, w_mt = [], []
                    for z in AA:
                        if len(AA[z]) >= 2:
                            change = []
                            for h in range(0,len(AA[z])):
                                gene_aj, rf_aj, at_aj, nt_aj, aa_aj, tps_aj = vsg.split('_')[0], vsg.split('_')[1], ALT[vsg.split('_')[1]], AA[z][h].split('\t')[0], AA[z][h].split('\t')[1], AA[z][h].split('\t')[2]
                                pos_aj = get_pos(vsg, tps_aj)
                                change.append(gene_aj+'\t'+pos_aj+'\t'+rf_aj+'\t'+at_aj+'\t'+nt_aj+'\t'+aa_aj+'\t'+tps_aj)
                            w_a_new, w_mt_new = get_merge(gene_aj, change)
                            w_a.append(w_a_new)
                            w_mt.append(w_mt_new)
                        else:
                            w_a.append(AA[z][0].split('\t')[1])
                            w_mt.append(AA[z][0].split('\t')[3])
                    if flg > 0:
                        var_out += var_inf[j]+'|'+':'.join(w_c)+'|'+':'.join(w_a)+'|'+':'.join(w_mt)+';'
                    else:
                        var_out += var_inf[j]+'|NA|NA|NA;'
            print var_out
                
def get_pos(sgf, wpsf):
    psf, rff, tpsf = sgf.split('_')[2], sgf.split('_')[1], sgf.split('|')[1]
    ps_nd = int(psf)+int(wpsf)-int(tpsf)
    if rff == 'T' or rff == 'G':
        int(psf)-int(wpsf)+int(tpsf)
    ps_nd = str(ps_nd)
    return ps_nd

def get_merge(ge, carray):
    mergent, mergeaa, mergetp, mergeot = '', '', '', ''
    NT = {'A':'T','T':'A','G':'C','C':'G'}
    strand, ch_mg = '+', '13'
    raa, raam = 'NA', 'NA'
    if ge == 'BRCA1':
        strand = '-'
        ch_mg = '17'
    pos_mg, nt_mg, at_mg, tps_mg, cc_mg, aa_mg = [], [], [], [], [], []
    for n in range(0, len(carray)):
        inf_car = carray[n].split('\t')
        pos_mg.append(int(inf_car[1]))
        if 'c.' in inf_car[4]:
            nt_mg.append(int(re.search(r'(\d+)',inf_car[4]).group()))
        else:
            nt_mg.append(0)
        cc_mg.append(inf_car[4])
        aa_mg.append(inf_car[5])
        tps_mg.append(int(inf_car[6]))
        if 'p.' in inf_car[5]:
            raa = re.search(r'([A-Z])\d+',inf_car[5]).group(1)
            raam = re.search(r'([A-Z])\d+',inf_car[5]).group(0)
        if strand == '+':
            at_mg.append(inf_car[3])
        else:
            at_mg.append(NT[inf_car[3]])
    tps_mg.sort()
    tps_mg = map(str, tps_mg)
    mergeot = ':'.join(tps_mg)+'|'+':'.join(cc_mg)+'|'
    if raa == 'NA':
        mergeot += aa_mg[0]+'|'+aa_mg[0]+';'
        mergent = 'NA'
        mergeaa, mergetp = merge_mut(mergent, raa)
    elif len(carray) < 3:
        mergent = merge_nt(nt_mg, ch_mg, pos_mg, at_mg, strand)
        mergeaa, mergetp = merge_mut(mergent, raa)
        mergeot += raam+mergeaa+'|'+mergetp+';'
    else:
        mergent = ''.join(at_mg)
        mergeaa, mergetp = merge_mut(mergent, raa)
        mergeot += raam+mergeaa+'|'+mergetp+';'
        split1, split2 = [0, 0, 1], [1, 2, 2]
        nt_mg_tmp, pos_mg_tmp, at_mg_tmp, tps_mg_tmp, cc_mg_tmp = [], [], [], [], []
        for l in  range(0, len(split1)):
            nt_mg_tmp = [nt_mg[split1[l]], nt_mg[split2[l]]]
            pos_mg_tmp = [pos_mg[split1[l]], pos_mg[split2[l]]]
            at_mg_tmp = [at_mg[split1[l]], at_mg[split2[l]]]
            tps_mg_tmp = [tps_mg[split1[l]], tps_mg[split2[l]]]
            cc_mg_tmp = [cc_mg[split1[l]], cc_mg[split2[l]]]
            mergent = merge_nt(nt_mg_tmp, ch_mg, pos_mg_tmp, at_mg_tmp, strand)
            mergeaa, mergetp = merge_mut(mergent, raa)
            mergeot += ':'.join(tps_mg_tmp)+'|'+':'.join(cc_mg_tmp)+'|'+raam+mergeaa+'|'+mergetp+';'

    return 'p.'+raam+mergeaa, mergetp

    
def merge_mut(mgnt, ra):
    AAchange = {'TTT' : 'F','TTC' : 'F','TTA' : 'L','TTG' : 'L','TCT' : 'S','TCC' : 'S','TCA' : 'S','TCG' : 'S','TAT' : 'Y','TAC' : 'Y','TAA' : 'X','TAG' : 'X','TGT' : 'C','TGC' : 'C','TGA' : 'X','TGG' : 'W','CTT' : 'L','CTC' : 'L','CTA' : 'L','CTG' : 'L','CCT' : 'P','CCC' : 'P','CCA' : 'P','CCG' : 'P','CAT' : 'H','CAC' : 'H','CAA' : 'Q','CAG' : 'Q','CGT' : 'R','CGC' : 'R','CGA' : 'R','CGG' : 'R','ATT' : 'I','ATC' : 'I','ATA' : 'I','ATG' : 'M','ACT' : 'T','ACC' : 'T','ACA' : 'T','ACG' : 'T','AAT' : 'N','AAC' : 'N','AAA' : 'K','AAG' : 'K','AGT' : 'S','AGC' : 'S','AGA' : 'R','AGG' : 'R','GTT' : 'V','GTC' : 'V','GTA' : 'V','GTG' : 'V','GCT' : 'A','GCC' : 'A','GCA' : 'A','GCG' : 'A','GAT' : 'D','GAC' : 'D','GAA' : 'E','GAG' : 'E','GGT' : 'G','GGC' : 'G','GGA' : 'G','GGG' : 'G', 'NA' : 'NA'}
    mgaa = AAchange[mgnt]
    mgtp = ''
    if mgaa == ra:
        mgtp = 'synonymous SNV'
    elif mgaa == 'X':
        mgtp = 'stopgain'
    elif mgaa == 'NA':
        mgtp = 'NA'
    else:
        mgtp = 'nonsynonymous SNV'
    return mgaa, mgtp

def merge_nt(ntarray, ch, posarray, atarray, strd):
    reference = '/biocluster/data/bioexec/database/genome/Homo_sapiens_assembly19/v1/Homo_sapiens_assembly19.fasta'
    NT_tmp = {'A':'T','T':'A','G':'C','C':'G'}
    mnt = ''
    nt_max = max(ntarray)
    ps_tmp, at_tmp = ch+':', ''
    if nt_max % 3 == 0:
        if abs(ntarray[0]-ntarray[1]) == 1:
            if strd == '+':
                ps_tmp = ps_tmp+str(posarray[ntarray.index(nt_max)]-2)+'-'+str(posarray[ntarray.index(nt_max)]-2)
                #print ps_tmp
                at_tmp = os.popen("samtools faidx  %s %s" % (reference, ps_tmp)).read().strip().split('\n')[1]
                mnt = at_tmp+atarray[0]+atarray[1]
            else:
                ps_tmp = ps_tmp+str(posarray[ntarray.index(nt_max)]+2)+'-'+str(posarray[ntarray.index(nt_max)]+2)
                at_tmp = NT_tmp[os.popen("samtools faidx  %s %s" % (reference, ps_tmp)).read().strip().split('\n')[1]]
                mnt = at_tmp+atarray[0]+atarray[1]
        else:
            ps_tmp = ps_tmp+str(int(np.mean(posarray)))+'-'+str(int(np.mean(posarray)))
            at_tmp = os.popen("samtools faidx  %s %s" % (reference, ps_tmp)).read().strip().split('\n')[1]
            if strd == '-':
                at_tmp = NT_tmp[at_tmp]
            mnt = atarray[0]+at_tmp+atarray[1]
    else:
        if strd == '+':
            ps_tmp = ps_tmp+str(posarray[ntarray.index(nt_max)]+1)+'-'+str(posarray[ntarray.index(nt_max)]+1)
            at_tmp = os.popen("samtools faidx  %s %s" % (reference, ps_tmp)).read().strip().split('\n')[1]
            mnt = atarray[0]+atarray[1]+at_tmp
        else:
            ps_tmp = ps_tmp+str(posarray[ntarray.index(nt_max)]-1)+'-'+str(posarray[ntarray.index(nt_max)]-1)
            at_tmp = NT_tmp[os.popen("samtools faidx  %s %s" % (reference, ps_tmp)).read().strip().split('\n')[1]]
            mnt = atarray[0]+atarray[1]+at_tmp
    #if mnt == '':
    #    print carray
    #    print str(nt_max)+'\t'+strand
    return mnt    



          
def read_anfile(anfile):
    An = {}
    with open(anfile) as af:
        for aline in af:
            ainf = aline.strip().split('\t')
            asg, anm, amtype, anno = ainf[0], ainf[1], ainf[2], ainf[3]
            anno_inf = anno.split(';')
            An.setdefault(asg, {})
            An[asg]['NM'] = anm
            An[asg]['mtype'] = amtype
            An[asg].setdefault('anno',[])
            for i in range(0, len(anno_inf)-1):
                An[asg]['anno'].append(anno_inf[i])

    return An



if __name__ == "__main__":
    main()
