count=$1
sample=$2
control=$3
efficiency=$4
lib=$5

bin='IdentifyLOFVariants'

#median normalized
python $bin/median_ratio.py -f $count > $sample\_sgRNA_count_medianratio.txt 
python $bin/median_ratio.py -f $control > control_sgRNA_count_medianratio.txt

#calculate growth rate by fitting    
python $bin/growth_ratio.py -f control_sgRNA_count.medianratio.txt -o control_sgRNA_count.median > control_sgRNA_count.growth.txt

#calculate score
python $bin/score.py -f $sample\_sgRNA_count_medianratio.txt -e $efficiency -g control_sgRNA_count.growth.txt -l $lib > $sample\.score.xls

#GMM fit
python $bin/gmm_fit.py -f $sample\.score.xls -n 4 -o $sample\_gmm > $sample\.out
python $bin/get_pvalue.py -f $sample\_gmm_pred.txt > $sample\_gmm_pred-addp.xls
