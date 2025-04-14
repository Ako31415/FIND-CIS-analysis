#!/bin/bash
#run by ./hallifted_bed_to_count.sh genotype dir Indel_pos_-start_stop_maternal_coodinate Indel_pos_-start_stop_paternal_coodinate genome_size_file bedgraph_condition1 bedgraph_condition2 cond1 cond2
EXPECTED_ARGS=9
E_BADARGS=8
if [ $# -ne $EXPECTED_ARGS ]
then
  echo "Usage: ./hallifted_bed_to_count.sh genotype dir Indel_pos_-start_stop_maternal_coodinate Indel_pos_-start_stop_paternal_coodinate   bedgraph_condition1(including 0 entries) bedgraph_condition2(including 0 entries) cond1 cond2 "
  exit $E_BADARGS
fi


g=$1
dir=$2
ID_M=$3
ID_P=$4
gen=$5
bed1=$6
bed2=$7
C1=$8
C2=$9


mkdir ${dir}
 cd ${dir}
 
 
#Indels were already cleared of duplicates by Amelie and also only those retained that are in two lines

#Since this script is for whole Indels and 3bp+- for no insertion, we need to make a bedgraph that is per base first to be able to count excatly. We only want this for the bases we need to save time and space.

echo "remove lines that contain an insertion larger than length indicated (i.e. a structural variant)"


gawk -v OFS='\t' '{split($4,var,"_");split(var[1],AL,":");if($3-$2==6){print $0}else{if($3-$2==length(AL[1])-1){print$0}else{if($3-$2==length(AL[2])-1){print$0}}}}' ${ID_P} > ${g}.${g}_coord.cleaned.bed


# echo "make bedfile containing regions that will be analysed"
# 
# cat ${ID_M} ${g}.${g}_coord.cleaned.bed| sortBed -g ${gen}| mergeBed > Merged_3bp.bed

echo "map counts on Indels and corresponding positions "
 
 zcat ${bed1} | intersectBed -a ${g}.${g}_coord.cleaned.bed -b - -wao | gawk -v OFS='\t' '{if(NR==1){CHR=$1;St=$2;Sp=$3;ID=$4;n=$9;AddVal=$8*$9}else{if($4==ID){n=n+$9;AddVal=AddVal+($8*$9)}else{print CHR,St,Sp,ID,AddVal/n;CHR=$1;St=$2;Sp=$3;ID=$4;n=$9;AddVal=$8*$9}}}END{print CHR,St,Sp,ID,AddVal/n}' > ${g}.B73.${C1}.counts.bed
 
 zcat ${bed1} | intersectBed -a ${ID_M} -b - -wao | gawk -v OFS='\t' '{if(NR==1){CHR=$1;St=$2;Sp=$3;ID=$4;n=$9;AddVal=$8*$9}else{if($4==ID){n=n+$9;AddVal=AddVal+($8*$9)}else{print CHR,St,Sp,ID,AddVal/n;CHR=$1;St=$2;Sp=$3;ID=$4;n=$9;AddVal=$8*$9}}}END{print CHR,St,Sp,ID,AddVal/n}' > B73.${g}.${C1}.counts.bed

 zcat ${bed2} | intersectBed -a ${g}.${g}_coord.cleaned.bed -b - -wao | gawk -v OFS='\t' '{if(NR==1){CHR=$1;St=$2;Sp=$3;ID=$4;n=$9;AddVal=$8*$9}else{if($4==ID){n=n+$9;AddVal=AddVal+($8*$9)}else{print CHR,St,Sp,ID,AddVal/n;CHR=$1;St=$2;Sp=$3;ID=$4;n=$9;AddVal=$8*$9}}}END{print CHR,St,Sp,ID,AddVal/n}' > ${g}.B73.${C2}.counts.bed
 
 zcat ${bed2} | intersectBed -a ${ID_M} -b - -wao | gawk -v OFS='\t' '{if(NR==1){CHR=$1;St=$2;Sp=$3;ID=$4;n=$9;AddVal=$8*$9}else{if($4==ID){n=n+$9;AddVal=AddVal+($8*$9)}else{print CHR,St,Sp,ID,AddVal/n;CHR=$1;St=$2;Sp=$3;ID=$4;n=$9;AddVal=$8*$9}}}END{print CHR,St,Sp,ID,AddVal/n}' > B73.${g}.${C2}.counts.bed

echo "cleaning up"

rm ${g}.${g}_coord.cleaned.bed


