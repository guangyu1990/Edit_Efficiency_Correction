step_anno=$1  #BRCA_NGG_sgRNA_step.format.txt
sample=$2
lib=$3 #CG_NGG/AT_NGG/CG_NG/AT_NG
contrl=$4
vardir=$outdir/group_var_call
db='Database'

bin='OutComeMap'

#co-edit
python $bin/co_edit.py -v $vardir/$sample\_var.out -s $db/BRCA_merg_oligo.fasta -m $db/BRCA_total_sgRNA_same.xls -o $sample\_codit.txt

#format window anno
python $bin/anno_window.py -mf $step_anno > $lib\_window_anno.txt

#anno lib sample
python $bin/co_anno_add.py -v $sample\_codit.txt -a $lib\_window_anno.txt > $sample\_codit_anno.txt

#pick outcome of lib sample
python $bin/map_window.py -v $sample\_codit_anno.txt > $sample\_codit_anno_pick.txt

#anno control sample
python $bin/co_anno_add.py -v $contrl -a $lib\_window_anno.txt > $lib\_CTRL_codit_anno.txt

#pick outcome of control sample
python $bin/map_window.py -v $lib\_CTRL_codit_anno.txt > $lib\_CTRL_codit_anno_pick.txt

#filt and format
python $bin/format_map.py -v $sample\_codit_anno_pick.txt -c $lib\_CTRL_codit_anno_pick.txt -o $sample\_codit_map.txt
