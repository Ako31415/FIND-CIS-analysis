#!/bin/bash
#run by ./counts_to_table.sh genotype dir cond1 cond2 Peaks_cond1 Peaks_cond2 bedgraph_condition1 bedgraph_condition2 genotype_file
EXPECTED_ARGS=9
E_BADARGS=9
if [ $# -ne $EXPECTED_ARGS ]
then
  echo "Usage: ./counts_to_table.sh genotype dir cond1 cond2 Peaks_cond1 Peaks_cond2 bedgraph_condition1 bedgraph_condition2 genotype_file "
  exit $E_BADARGS
fi


g=$1
dir=$2
C1=$3
C2=$4
PeakC1=$5
PeakC2=$6
bed1=$7
bed2=$8
geno=$9

cd ${dir}

#Start with : Chr start stop ID Value

#To add peak information

echo "sort and name peaks"

#Peaks are taken from Narrowpeak files of MACS3

for tr in ${C1};  do sort -k1,1 -k2,2n ${PeakC1} | gawk -v OFS='\t' '{if($2<$3 && $1!~"scaf-alt"){if($4~"peak"){print $1,$2,$3,$4"_"$5}else{print $1,$2,$3,"peak_"NR}}}' > ${g}.${tr}.peaks.names.sorted.bed ; done

for tr in ${C2};  do sort -k1,1 -k2,2n ${PeakC2} | gawk -v OFS='\t' '{if($2<$3 && $1!~"scaf-alt"){if($4~"peak"){print $1,$2,$3,$4"_"$5}else{print $1,$2,$3,"peak_"NR}}}' > ${g}.${tr}.peaks.names.sorted.bed ; done

echo "sort count files, mark heterozygous sites (scaf-alt)"

for tr in ${C1} ${C2}; do gawk -v OFS='\t' -v g=$g '{if($4!~"scaf-alt"){print $0} else{gsub("scaf-alt","scaf",$4);print $1,$2,$3,$4,"het"}}' ${dir}/B73.${g}.${tr}.counts.bed | sort -k1,1 -k2,2n  > B73.${g}.${tr}.counts.sorted.bed; done

for tr in ${C1} ${C2};  do gawk -v OFS='\t' -v g=$g '{if($1!~"scaf-alt"){print $0} else{gsub("scaf-alt","scaf",$1);print $1,$2,$3,$4,"het"}}' ${dir}/${g}.B73.${tr}.counts.bed | sort -k1,1 -k2,2n  > ${g}.B73.${tr}.counts.sorted.bed; done



#What we want for each genotype here: Chr  Start Stop   ID(INDEL,no,B73_pos)     Count  Peak_overlap   
echo "intersect count files with peaks"

for tr in ${C1} ${C2}; do intersectBed -a B73.${g}.${tr}.counts.sorted.bed -b ${g}.${tr}.peaks.names.sorted.bed -wa -wb -sorted -loj| sort  >  B73.${g}.${tr}.ID.count.sort.csv; done

for tr in ${C1} ${C2}; do intersectBed -a ${g}.B73.${tr}.counts.sorted.bed -b ${g}.${tr}.peaks.names.sorted.bed -wa -wb -sorted -loj| sort -k1,1 -k2,2n  > ${g}.B73.${tr}.ID.count.sort.csv; done


#Now make new file with both ref and alt, one WW, one DS, keep also ID, add genotype:

# What we want: Chr(B73)  STop=POS(B73 from ID)    REFal(B73)   ALTal(Paternal) ID_new(NAMchr_NAMPos)   Genotype (0/0 if B73, 1/1 if Pat) norm_Count(B73)   norm_Count(Pat) Peak(B73) Peak(Pat)  BF(BF (B73/(B73+Pat)), if both 0 print n.r.)  PEAK(B73) Peak(Pat)  

# the positions in this case should always be the last common base before the INDEL. This can be taken from the ID for B73 but needs to be calculated for NAM depending on if the line had the INDEL or not. Thus first make a file of NAM that contains genotype info and a NAM ID:

echo "adding genotype"

for tr in ${C1} ${C2}; do gawk -v OFS='\t' '{split($4, a, ":|_"); if( ((length(a[2])==1) && (a[7]=="del")) || ((length(a[2])>1) && (a[7]=="ins")) ){print substr($4,1,length($4)-4),a[7]}}' ${geno} | gawk -v OFS='\t' '{if(NR==FNR) {a[$1]=$1;t[$1]=$2; next} if($4 in a){if(t[$4]=="del"){print $0,"1/1",$1":"$2+3};if(t[$4]=="ins"){print $0,"1/1",$1":"$2}}else{split($4, var, ":|_");if(length(var[1])>1){print $0,"0/0",$1":"$2}else{print $0,"0/0",$1":"$2+3}}}' - ${g}.B73.${tr}.ID.count.sort.csv > ${g}.B73.${tr}.ID.count.sort.geno.csv;done


#now merge B73 and NAM files and calculate binding frequencies


echo "preparing binding frequency files"

for tr in ${C1} ${C2}; do gawk -v OFS='\t' '{if(NR==FNR) {a[$4]=$4;Value[$4]=$5;Peak[$4]=$9;geno[$4]=$10;ID[$4]=$11; next} if($4 in a){split($4,var,"_");split(var[1],AL,":");if(length(AL[1])>1){print $1,$2,AL[1],AL[2],ID[$4],geno[$4],$5,Value[$4],$9,Peak[$4]}else{print $1,$2+3,AL[1],AL[2],ID[$4],geno[$4],$5,Value[$4],$9,Peak[$4]}}}' ${g}.B73.${tr}.ID.count.sort.geno.csv B73.${g}.${tr}.ID.count.sort.csv |  gawk -v OFS='\t' '{if($7+$8>0){if($9~"peak" || $10~"peak"){print $0,($7/($7+$8))}else{print $0,"n.p."($7/($7+$8))}}else{print $0,"n.p.n.r."}}'  > B73.${g}.${tr}.PF.GT.csv; done


echo "adding read depth information "


#Chr(B73)  STop=POS(B73) REFal(B73)  ALTal(Pat) ID(Patchr_PatPOS) Genotype (0/0 if B73, 1/1 if Pat) norm_Count(B73)   norm_Count(Pat) PEAK(B73) Peak(Pat)  BF (n.p. if no peak in either B73 or Pat) read_Count(B73) read_Count(Pat)

#Determine value per read first

CW=`zcat ${bed1} | awk -v OFS="\t" 'BEGIN{n=10000} {if($4>0 && $4<n){n=$4}} END { print n }' `;
CD=`zcat ${bed2} | awk -v OFS="\t" 'BEGIN{n=10000} {if($4>0 && $4<n){n=$4}} END { print n }' `;

for tr in ${C1} ; do gawk -v OFS='\t' -v g=$g -v CW=${CW} '{{print $0,int(($7/CW)+0.5),int(($8/CW)+0.5)}}' B73.${g}.${tr}.PF.GT.csv  > B73.${g}.${tr}.PF.GT.RN.csv; done

for tr in ${C2} ; do gawk -v OFS='\t' -v g=$g -v CD=${CD} '{{print $0,int(($7/CD)+0.5),int(($8/CD)+0.5)}}' B73.${g}.${tr}.PF.GT.csv  > B73.${g}.${tr}.PF.GT.RN.csv; done


rm B73.${g}.${C1}.counts.sorted.bed

rm ${g}.B73.${C1}.counts.sorted.bed


rm B73.${g}.${C1}.ID.count.sort.csv

rm B73.${g}.${C2}.ID.count.sort.csv

rm ${g}.B73.${C1}.ID.count.sort.csv

rm ${g}.B73.${C2}.ID.count.sort.csv

rm B73.${g}.${tr}.PF.GT.csv
