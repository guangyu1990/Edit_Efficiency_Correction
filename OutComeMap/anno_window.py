#!/usr/bin/python

import click
import re
import os
import math
import pandas as pd
import numpy as np
import scipy
from scipy import stats

@click.command()
@click.option("--mutfile","-mf",help="mut anno file")

def main(mutfile):
    Window = {}
    position = 'BRCA\t0'
    with open(mutfile) as mf:
        for mline in mf:
            if mline.startswith('Chr'):
                continue
            minf = mline.strip().split('\t')
            msg, mpos, mrw_nt, mrw_aa = minf[0], minf[2], minf[8], minf[9]
            if mrw_nt == '.':
                minf[8] = 'NA'
            if mrw_aa == '.':
                minf[9] = 'NA'
            mline_new = '\t'.join(minf[0:8])+'\t'+minf[8]+'\t'+minf[9]+'\t'+'\t'.join(minf[10:len(minf)])
            if msg == position.split('\t'[0]) and abs(int(mpos)-int(position.split('\t')[1])) <= 4:
                Window[msg].append(mline_new)
            else:
                Window.setdefault(msg, [])
                Window[msg].append(mline_new)
                position = msg+'\t'+mpos

    for ms in Window:
        if len(Window[ms]) >= 2:
            #print ms
            AA = {}
            MAA = []
            maa = ''
            out_w = ''
            w_anno = ''
            for i in range(0,len(Window[ms])):
                inf_w = Window[ms][i].split('\t')
                if inf_w[0].split('_')[2] == inf_w[2]:
                    out_w = inf_w      
                aa = inf_w[9]          
                if aa == 'NA':
                    aa = inf_w[7]
                    MAA.append(aa)
                    AA.setdefault(aa,[])
                    AA[aa].append(Window[ms][i])
                    w_anno += inf_w[12]+'|'+inf_w[8]+'|'+inf_w[7]+'|'+inf_w[7]+';'
                    continue
                w_anno += inf_w[12]+'|'+inf_w[8]+'|'+inf_w[9]+'|'+inf_w[7]+';'
                aa = re.search(r'([A-Z]\d+)',aa).group()
                AA.setdefault(aa,[])
                AA[aa].append(Window[ms][i])
                #print Window[ms][i]
            for a in AA:
                if len(AA[a]) >= 2:
                    change = []
                    for j in range(0,len(AA[a])):
                        inf_aj = AA[a][j].split('\t')
                        ch_aj, gene_aj, type_aj = inf_aj[1], inf_aj[5], inf_aj[7]
                        pos_aj, rf_aj, at_aj, nt_aj, aa_aj, tps_aj = inf_aj[2], inf_aj[3], inf_aj[4], inf_aj[8], inf_aj[9], inf_aj[12]
                        if aa_aj == 'NA':
                            aa_aj = type_aj
                        change.append(ch_aj+'\t'+pos_aj+'\t'+rf_aj+'\t'+at_aj+'\t'+nt_aj+'\t'+aa_aj+'\t'+tps_aj)
                        MAA.append(type_aj)
                    maa_merge, anno_merge = get_merge(gene_aj, change)
                    w_anno += anno_merge
                    MAA.append(maa_merge)
                else:
                    inf_a = AA[a][0].split('\t')
                    MAA.append(inf_a[7])
            if 'stopgain' in MAA:
                maa = 'stopgain'
            elif 'stoploss' in MAA:
                maa = 'stoploss'
            elif 'splicing' in MAA:
                maa = 'splicing'
            elif 'nonsynonymous SNV' in MAA:
                maa = 'nonsynonymous SNV'
            elif 'UTR5' in MAA:
                maa = 'UTR5'
            elif 'UTR3' in MAA:
                maa = 'UTR3'
            elif 'synonymous SNV' in MAA:
                maa = 'synonymous SNV'
            else:
                maa = MAA[0]
            print out_w[0]+'\t'+out_w[6]+'\t'+maa+'\t'+w_anno
        else:
            out_t = Window[ms][0].split('\t')
            if out_t[9] == 'NA':
                out_t[9] = out_t[7]
            print out_t[0]+'\t'+out_t[6]+'\t'+out_t[7]+'\t'+out_t[12]+'|'+out_t[8]+'|'+out_t[9]+'|'+out_t[7]+';'
            #maa = ''

def get_merge(ge, carray):
    mergent, mergeaa, mergetp, mergeot = '', '', '', ''
    reference = '/public/database/hg19/ucsc.hg19.fasta'
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
        mergeot += 'p.'+raam+mergeaa+'|'+mergetp+';'
    else:
        mergent = ''.join(at_mg)
        mergeaa, mergetp = merge_mut(mergent, raa)
        mergeot += 'p.'+raam+mergeaa+'|'+mergetp+';'
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
            mergeot += ':'.join(tps_mg_tmp)+'|'+':'.join(cc_mg_tmp)+'|'+'p.'+raam+mergeaa+'|'+mergetp+';'

    return mergetp, mergeot

    
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

if __name__ == "__main__":
    main()
