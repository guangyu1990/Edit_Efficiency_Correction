#!/usr/bin/python
import re
import click

@click.command()
@click.option("--vrfile","-v",help="var file")
@click.option("--cnfile","-c",help="control var file", default='NA')
@click.option("--lib","-l",help="lib type")
@click.option("--outfile","-o",help="outfile")

def main(vrfile, cnfile, lib, outfile):
    NT = {'A':'G', 'T':'C', 'G':'A', 'C':'T'}
    CN, Target = {}, {}
    if cnfile != 'NA':
        CN = read_cnfile(cnfile)
    Var = read_vrfile(vrfile)
    var_out = 'sgRNA\tpos\tref\talt\tref_alt\tvars\tfreq\tflag\n'
    open(outfile,'w').writelines(var_out)
    OT = open(outfile,'a+')
    for vs in Var:
        vs_inf = vs.split('\t')
        sg, ref, pos, dep = vs_inf[0], vs_inf[1], vs_inf[2], int(vs_inf[3])
        dep_tmp, var, freq = dep, 'NA', 0 
        slib, tpos = '', sg.split('|')[1]
        if '_NG_' in sg:
            slib = 'NG'
        else:
            slib = 'NGG'
        if ref == 'C':
            slib = 'CG_'+slib
        elif ref == 'A':
            slib = 'AT_'+slib
        #print slib
        if slib == lib and int(tpos)+20 == int(pos):
            for va in Var[vs]:
                freq_tmp = Var[vs][va].split('\t')[0]
                dep = int(dep - dep_tmp*float(freq_tmp))
                if NT[ref] == va:
                    var = va
                    freq = dep_tmp*float(freq_tmp)
        if var != 'NA':
            error = 0
            if sg in CN:
                #print sg
                error = (freq+dep)*float(CN[sg])
            var_ref = int(freq+dep-error)
            freq = freq-error
            if freq <= 0 or var_ref == 0:
                freq = 0
            else:
                freq = freq/var_ref
            var_out = sg+'\t'+pos+'\t'+ref+'\t'+var+'\t'+str(var_ref)+'\t'+str(freq)+'\ttarget\n' 
            OT.writelines(var_out)
    OT.close()
            

def read_vrfile(vfile):
    V = {}
    with open(vfile) as vf:
        for vline in vf:
            if vline.startswith('sgRNA'):
                continue
            vinf = vline.split('\t')
            vsg, vpos, vref, valt, vdep, vfreq, vflag = vinf[0], vinf[1], vinf[2], vinf[3], vinf[4], vinf[5], vinf[6]
            vid = vsg+'\t'+vref+'\t'+vpos+'\t'+vdep
            V.setdefault(vid, {})
            V[vid][valt] = vfreq+'\t'+vflag
    return V

def read_cnfile(cfile):
    C = {}
    with open(cfile) as cf:
        for cline in cf:
            if cline.startswith('sgRNA'):
                continue
            cinf = cline.strip().split('\t')
            csg, cdp, cfq = cinf[0], cinf[4], cinf[5]
            if float(cdp) >= 50:
                C[csg] = cfq
    return C

if __name__=='__main__':
   main()
