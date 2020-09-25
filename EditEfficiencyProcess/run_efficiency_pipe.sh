bam=$1
sample=$2
lib=$3
outdir=$4
groupdir=$outdir/read_group
vardir=$outdir/group_var_call
bin='EditEfficiencyProcess'
db='Database'

#bam2fasta
python $bin/bam_subsgRNA.py -b $bam -o $groupdir/$sample >  $groupdir/$sample\_bam.out

#flt recombine reads
python $bin/recombine_flt.py -q1 $groupdir/$sample\_sgseq-first.fasta -q2 $groupdir/$sample\_sgseq-second.fasta -i 0.75 -d jaro -o $groupdir/$sample > $groupdir/$sample\_fltrb.out

#group to ref sgRNA seq
cat $groupdir/$sample\_sgseq-first_fltrb.fasta $db/BRCA_merg_oligo.fasta > $groupdir/$sample\_sgseq-first_addseq.fasta
cd-hit -i $groupdir/$sample\_sgseq-first_addseq.fasta -o $groupdir/$sample\_cdhit_85.out -c 0.85 -n 5 -M 30000 -T 16 -d 0

#group2var-read
python $bin/group2varead.py -g $groupdir/$sample\_cdhit_85.out.clstr -s $db/BRCA_merg_oligo-uniq.txt -q1 $groupdir/$sample\_sgseq-first_fltrb.fasta -q2 $groupdir/$sample\_sgseq-second_fltrb.fasta -l $lib -o $vardir/$sample\_var_read.txt > $vardir/$sample\_read.out

#var call
python $bin/varead2var.py -v $vardir/$sample\_var_read.txt -s $db/BRCA_total_sgRNA_same.xls -l $lib -o $vardir/$sample\_var.txt > $vardir/$sample\_var.out 
