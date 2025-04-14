#!/bin/bash
#run by ./get_control_counts.sh genotype dir window cond1 cond2 SNPs_tested_bino_cond1 genome_size_file sig_SNPs_cond2 control_bedgraph
EXPECTED_ARGS=9
E_BADARGS=9
if [ $# -ne $EXPECTED_ARGS ]
then
  echo "Usage: ./get_control_counts.sh genotype dir window(even number) cond1 cond2 SNPs_tested_bino_cond1 genome_size_file sig_SNPs_cond2 control_bedgraph "
  exit $E_BADARGS
fi


g=$1
dir=$2
window=$3
cond1=$4
cond2=$5
bino=$6
genome=$7
bino2=$8
bedG=$9



cd ${dir}

#Start with : Chr start stop ID Value

#To add peak information

echo "get ${window}bp around INDEL ${cond1}"

#Since the base that is recorded is the last common one, this base minus half the window is the start site. 
gawk -v OFS='\t' -v g=$g -v window=$window '{if(NR>1){print substr($3,2,length($3)-2),$4-(int(window/2)),$4+(int(window/2)-1)+(length($5)-2),substr($7,2,length($7)-2)}}' ${bino} |  gawk -v OFS='\t' '{if($2<0){print $1,0,$3,$4}else{print $0}}' | sortBed -g ${genome} > B73.${g}.${cond1}.${window}.bed


gawk -v OFS='\t' -v g=$g -v window=$window '{if(NR>1){split($7,pos,":"); print substr(pos[1],2,length(pos[1])-1),substr(pos[2],1,length(pos[2])-1)-(int(window/2)),substr(pos[2],1,length(pos[2])-1)+(int(window/2)-1)+(length($6)-2),substr($7,2,length($7)-2)}}' ${bino}  |  gawk -v OFS='\t' '{if($2<0){print $1,0,$3,$4}else{print $0}}'| sortBed -g ${genome} > ${g}.${cond1}.${window}.bed

echo "get ${window}bp around SNP ${cond2}"

gawk -v OFS='\t' -v g=$g -v window=$window '{if(NR>1){print substr($3,2,length($3)-2),$4-(int(window/2)),$4+(int(window/2)-1)+(length($5)-2),substr($7,2,length($7)-2)}}' ${bino2} |  gawk -v OFS='\t' '{if($2<0){print $1,0,$3,$4}else{print $0}}' | sortBed -g ${genome} > B73.${g}.${cond2}.${window}.bed

 gawk -v OFS='\t' -v g=$g -v window=$window '{if(NR>1){split($7,pos,":"); print substr(pos[1],2,length(pos[1])-1),substr(pos[2],1,length(pos[2])-1)-(int(window/2)),substr(pos[2],1,length(pos[2])-1)+(int(window/2)-1)+(length($6)-2),substr($7,2,length($7)-2)}}' ${bino2}  |  gawk -v OFS='\t' '{if($2<0){print $1,0,$3,$4}else{print $0}}'| sortBed -g ${genome} > ${g}.${cond2}.${window}.bed
 
 
echo "get merged file both conditions"

cat ${g}.${cond1}.${window}.bed ${g}.${cond2}.${window}.bed | sortBed -g ${genome}| mergeBed > ${g}.both.${window}.bed

cat B73.${g}.${cond1}.${window}.bed B73.${g}.${cond2}.${window}.bed | sortBed -g ${genome}| mergeBed > B73.${g}.both.${window}.bed

 
 
 echo "get part of begraph that overlaps with features in question paternal allele"
 
 zcat ${bedG} | intersectBed -a - -b ${g}.both.${window}.bed -wa | sortBed -g ${genome}|uniq  > ${g}.both.bedgraph

echo "get part of begraph that overlaps with features in question maternal allele"
 
 zcat ${bedG} | intersectBed -a - -b B73.${g}.both.${window}.bed -wa | sortBed -g ${genome}|uniq  > B73.${g}.both.bedgraph

 
 
echo "get counts paternal allele"

gawk -v OFS='\t' '{{a=($3 - $2); if(a>1){for(i=a; i>=1; i--){ print $1, ($3-i),($3-i+1), $4}} else{ print $0 }}}' ${g}.both.bedgraph| bedtools map -a ${g}.${cond2}.${window}.bed -b - -g ${genome} -c 4 -o mean > ${g}.${cond2}.CTRL.${window}.bed


gawk -v OFS='\t' '{{a=($3 - $2); if(a>1){for(i=a; i>=1; i--){ print $1, ($3-i),($3-i+1), $4}} else{ print $0 }}}' ${g}.both.bedgraph| bedtools map -a ${g}.${cond1}.${window}.bed -b - -g ${genome} -c 4 -o mean > ${g}.${cond1}.CTRL.${window}.bed

echo "get counts maternal allele"

gawk -v OFS='\t' '{{a=($3 - $2); if(a>1){for(i=a; i>=1; i--){ print $1, ($3-i),($3-i+1), $4}} else{ print $0 }}}' B73.${g}.both.bedgraph| bedtools map -a B73.${g}.${cond2}.${window}.bed -b - -g ${genome} -c 4 -o mean > B73.${g}.${cond2}.CTRL.${window}.bed


gawk -v OFS='\t' '{{a=($3 - $2); if(a>1){for(i=a; i>=1; i--){ print $1, ($3-i),($3-i+1), $4}} else{ print $0 }}}' B73.${g}.both.bedgraph| bedtools map -a B73.${g}.${cond1}.${window}.bed -b - -g ${genome} -c 4 -o mean > B73.${g}.${cond1}.CTRL.${window}.bed

