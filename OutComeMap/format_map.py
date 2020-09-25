#!/usr/bin/python
import re
import click

@click.command()
@click.option("--vrfile","-v",help="var file")
@click.option("--cnfile","-c",help="control var file", default='NA')
@click.option("--outfile","-o",help="outfile")

def main(vrfile, cnfile, outfile):
    NT = {'A':'G', 'T':'C', 'G':'A', 'C':'T'}
    CN = {}
    if cnfile != 'NA':
        CN = read_cnfile(cnfile)
    Var = read_vrfile(vrfile)
    var_out = 'sgRNA\tpos\tref\talt\tref_alt\tvars\tfreq\taachange\tmutype\n'
    open(outfile,'w').writelines(var_out)
    OT = open(outfile,'a+')
    for vs in Var:
        vs_inf = Var[vs].split('\t')
        sg, pos, dep, var = vs, [], float(vs_inf[3]), float(vs_inf[4])
        ref, alt, aa, mtype = vs_inf[1], vs_inf[2], vs_inf[5], vs_inf[6]
        if '|' in vs_inf[0]:
            pos = vs_inf[0].split('|')
        else:
            pos.append(vs_inf[0])
        error = 0
        if sg in CN:
            cn_dep, cn_inf = float(CN[sg].split('\t')[0]), CN[sg].split('\t')[1].split(';')
            cn_var = 0
            for i in range(0, len(pos)):
                for j in range(0, len(cn_inf)-1):
                    cn_pos = cn_inf[j].split('|')[0]
                    if cn_pos == pos[i]:
                        cn_var += float(cn_inf[j].split('|')[1])
            error = float(cn_var)/float(cn_dep)
        dep = float(dep - dep*error)
        var = float(var - var*error)
        freq = 'NA'
        if dep <= 0 or var < 0:
            dep, var, freq = 0, 0, 0
        else:
            freq = float(var)/float(dep)
        var_out = sg+'\t'+'|'.join(pos)+'\t'+ref+'\t'+alt+'\t'+str(dep)+'\t'+str(freq)+'\t'+aa+'\t'+mtype+'\n' 
        OT.writelines(var_out)
    OT.close()
            

def read_vrfile(vfile):
    V = {}
    NT1 = {'A':'A', 'T':'A', 'C':'C', 'G':'C'}
    NT2 = {'A':'G', 'C':'T'}
    with open(vfile) as vf:
        for vline in vf:
            if vline.startswith('sgRNA'):
                continue
            vinf = vline.split('\t')
            vsg, vpos, vref, valt, vdep, vvar, vaa, vtyp = vinf[0], vinf[4], NT1[vinf[0].split('_')[1]], NT2[NT1[vinf[0].split('_')[1]]], vinf[2], vinf[7], vinf[5], vinf[6]
            vid = vsg
            V[vid] = vpos+'\t'+vref+'\t'+valt+'\t'+vdep+'\t'+vvar+'\t'+vaa+'\t'+vtyp
    return V

def read_cnfile(cfile):
    C = {}
    with open(cfile) as cf:
        for cline in cf:
            if cline.startswith('sgRNA'):
                continue
            cinf = cline.strip().split('\t')
            csg, cdep, cvar = cinf[0], cinf[2], cinf[11]
            if float(cdep) >= 50:
                C[csg] = cdep+'\t'+cvar
    return C

if __name__=='__main__':
   main()
