# Analysis pipeline for allele-specific MOA-seq analysis using INDELs instead of SNPs

*Main code was written by Amelie Kok. Scripts for integration with MOA-seq were written by Julia Engelhorn*

This part of the pipeline shows how a list of biallelic INDELs is created by pairwise genome alignments.

After mapping of MOA-seq reads to the genomes an analysis of allele-specific transcription factor (TF) binding can be performed similarly as for SNPs (described here: https://github.com/jengelhorn/AS-MOA ).

Additionally, an analysis of allele-specific methylation can be performed (see ../methylationAnalysis) and sites can also be used for bQTL mapping (see ../bQTL-Analysis).


Please note that INDEL detection by whole genome alignment is challenging in regions with e.g. local repeats or stretches of the same base. We noticed that especially positioning of one base pair INDELs in these regions can be ambiguous and lead to e.g. omission of some of these small deletions and hence an errouness genotype call. We included several quality control steps to exclude as many of these cases as possible, e.g. removal of overlapping INDLEs, translating the start and stop coordinates separately and re-checking the correct length of the insertion site. However, we cannot completely rule out that some ambiguous genotype information was retained in our final lists and therefore advise users to carefully check the local sequence when working with specific bQTL of interest.

<br/><br/>

## Step 1: pairwise genome alignments

Genomes of each NAM parent are aligned to one reference genome (here: B73) using Anchorwave v1.2.2.


<br/><br/>

1. Gene sequences of B73 are extracted using anchorwave gff2seq

```{bash}
anchorwave gff2seq -r Zm-B73-REFERENCE-NAM-5.0.B73.fa -i Zea_mays.Zm-B73-REFERENCE-NAM-5.0.57.B73-chr.gff3 -o B73_anchors_cds.fa
```

<br/><br/>

2. Gene sequences are aligned to the B73/reference sequence, and each of the other NAM parents genomes to localise anchors

```{bash}
minimap2 -x splice -t 50 -k 12 -a -p 0.4 -N 20 Zm-B73-REFERENCE-NAM-5.0.B73.fa B73_anchors_cds.fa > B73cds_against_B73.sam

for g in [NAMparents...]; do echo "${g}"; mkdir ${g}; minimap2 -x splice -t 10 -k 12 -a -p 0.4 -N 20 ${g}.pseudomolecules-v1-sm.fasta B73_anchors_cds.fa > B73cds_against_${g}.sam; done
```

<br/><br/>

3. Pairwise alignments are created using anchorwave proali

```{bash}
for g in [NAMparents...]; do echo "${g}"; anchorwave proali -i Zea_mays.Zm-B73-REFERENCE-NAM-5.0.57.B73-chr.gff3 -as B73_anchors_cds.fa -r Zm-B73-REFERENCE-NAM-5.0.B73.fa -a B73cds_against_${g}.sam -ar B73cds_against_B73.sam -s ${g}.pseudomolecules-v2.1-sm.fasta -n anchors -R 1 -Q 1 -o ${g}_againstB73.maf -t 14; done
```

<br/><br/>

## Step 2: Variant calling

Variants are extracted from the MAF (multiple alignment file) output files using wgatools wgatools 0.1.0.

```{bash}
for g in [NAMparents...]; do wgatools call ${g}_againstB73.maf -n ${g} -s -l 1 -t 80 -v > ${g}_againstB73.vcf; done
```

<br/><br/>

Chain files are created to enable lift over of coordinates between genomes:

```{bash}
for g in [NAMparents...]; do wgatools maf2chain ./${g}/${g}_againstB73.maf -t 55 -v > ./${g}/${g}_againstB73.chain; done
```

<br/><br/>

## Step 3: Determining biallelic INDELs

Create a BED file with only INDELs from the VCF files created in the previous step. Include the base before and after the indel, so as to have two bases on the deletion allele to serve as point of comparison.

```{bash}
for g in [NAMparents...]; do gawk -v OFS='\t' '{if(!($0 ~ "^#")){if((length($4)>1 && length($4)<51) || (length($5)>1 && length($5)<51)){print $1, $2-1, $2+length($4), $1":"$2":"$4":"$5}}}' ${g}_againstB73.vcf > ${g}_againstB73.INDELs.bed; done
```

Biallelic sites were determined as sites where every inbred line with the insertion allele had the same insertion sequence and position. The other allele remaining then was a shared deletion allele.

<br/><br/>

Create one file with the sequences of INDELs as 4th and 5th column:

```{bash}
for g in [NAMparents...]; do gawk -v OFS='\t' '{split($4, a, ":"); print $1, $2, $3, a[3]":"a[4]}' ${g}_againstB73.INDELs.bed >> All_againstB73.INDELs.seqs.tsv; done
```

<br/><br/>

Use bash sort and uniq -c to find out how often an INDEL occurs between B73 and each of the different NAM inbred lines. If an INDEL varies in sequence at a position there will be two separate entries.

```{bash}
sort All_againstB73.INDELs.seqs.tsv | uniq -c | gawk -v OFS='\t' '{print $1, $2, $3, $4, $5}' > All_againstB73.INDELs.seqs.uniq.tsv
```

<br/><br/>

Remove the count column in order to run bedtools intersect:

```{bash}
gawk -v OFS='\t' '{print $2, $3, $4, $5}' All_againstB73.INDELs.seqs.uniq.tsv | sort -k1,1 -k2,2n > All_againstB73.INDELs.seqs.uniq.noCounts.bed
```

<br/><br/>

Run bedtools intersect in order to identify INDELs overlapping in their B73 coordinates to exclude them:

```{bash}
bedtools intersect -a All_againstB73.INDELs.seqs.uniq.noCounts.bed -b All_againstB73.INDELs.seqs.uniq.noCounts.bed -wa -wb -loj > All_againstB73.INDELs.seqs.uniq.noCounts.intersect.tsv
```

<br/><br/>

Only keep entries that occur only once as overlapping with themselves and remove the rest: 

```{bash}
gawk -v OFS='\t' 'BEGIN{vars["-"]="-"; dups["-"]="-"}{if($1":"$2":"$3":"$4 in vars){ dups[$1":"$2":"$3":"$4]=$1":"$2":"$3":"$4 } else if($1":"$2":"$3":"$4==$5":"$6":"$7":"$8 || $5=="."){ vars[$1":"$2":"$3":"$4]=$1":"$2":"$3":"$4 } else { vars[$1":"$2":"$3":"$4]=$1":"$2":"$3":"$4 }} END{ for(i in vars){ if(!(i in dups)){ split(vars[i], a, ":"); print a[1], a[2], a[3], a[4]":"a[5] }}}' All_againstB73.INDELs.seqs.uniq.noCounts.intersect.tsv | sort -k1,1 -k2,2n > All_againstB73.INDELs.seqs.biallelic.bed
```

<br/><br/>

Take the counts of how often the INDELs occur in the different hybrids and add them to the biallic file in order to filter for a minimum of 2 hybrids the variants occur in:

```{bash}
gawk -v OFS='\t' '{if(NR==FNR){ counts[$2"_"$3"_"$4"_"$5]=$1; next} if($1"_"$2"_"$3"_"$4 in counts){ print $0, counts[$1"_"$2"_"$3"_"$4] }}' All_againstB73.INDELs.seqs.uniq.bed All_againstB73.INDELs.seqs.biallelic.bed | sort -k1,1 -k2,2n > All_againstB73.INDELs.seqs.biallelic.counts.tsv
```

<br/><br/>

Filter out positions that only occur in 1 hybrid:

```{bash}
gawk -v OFS='\t' '{if($5>1){ print $1, $2, $3, $4"_"$5 }}' All_againstB73.INDELs.seqs.biallelic.counts.tsv > All_againstB73.INDELs.seqs.biallelic.counts.min2hyb.bed
```

<br/><br/>

Create BED files with three positions before and after the INDEL. This is so that there are some bases to compare coverage for the deletion allele. Start and Stop coordinates are split into 4 separate rows to make liftover easier:

```{bash}
gawk -v OFS='\t' '{print $1, $2-2, $2-1, $4"_"$1"_"$2"_"$3"_St1\n" $1, $2, $2+1, $4"_"$1"_"$2"_"$3"_St2\n" $1, $3-1, $3, $4"_"$1"_"$2"_"$3"_Sp1\n" $1, $3+1, $3+2, $4"_"$1"_"$2"_"$3"_Sp2"}' All_againstB73.INDELs.seqs.biallelic.counts.min2hyb.bed > All_againstB73.INDELs.seqs.biallelic.counts.min2hyb.pm3bp.bed
```

<br/><br/>

Lift over files with start and stop positions of the selected INDELs to each NAM parent genome using CrossMap (version 0.7.0). The chain files created in step 2 are used.

```{bash}
for g in [NAMparents...]; do CrossMap bed --chromid a --unmap-file ${g}_againstB73.seqs.biallelic.pm3bp.unmapped.bed ${g}_againstB73.chain All_againstB73.INDELs.seqs.biallelic.counts.min2hyb.pm3bp.bed ${g}_againstB73.seqs.biallelic.pm3bp.crossmapped.bed; done
```

<br/><br/>

Merge start and stop positions accordingly after lift over, if all 4 coordinates are present:

```{bash}
for g in [NAMparents...]; do echo "${g}"; gawk -v OFS='\t' '{split($4, a, "_"); chr[a[1]"_"a[2]"_"a[3]"_"a[4]"_"a[5]]=$1; if(a[6]=="St1"){ st1[a[1]"_"a[2]"_"a[3]"_"a[4]"_"a[5]]=$3 } else if(a[6]=="St2"){ st2[a[1]"_"a[2]"_"a[3]"_"a[4]"_"a[5]]=$3 } else if(a[6]=="Sp1"){ sp1[a[1]"_"a[2]"_"a[3]"_"a[4]"_"a[5]]=$3 } else if(a[6]=="Sp2"){ sp2[a[1]"_"a[2]"_"a[3]"_"a[4]"_"a[5]]=$3 }} END { for(i in st1){ if((i in st2) && (i in sp1) && (i in sp2)) { if(st1[i]<st2[i]){ print chr[i], st1[i], st2[i], sp1[i], sp2[i], i } else if(st1[i]>st2[i]){ print chr[i], sp2[i], sp1[i], st2[i], st1[i], i }}}}' ${g}_againstB73.seqs.biallelic.pm3bp.crossmapped.bed | sort -k1,1 -k2,2n > ${g}_againstB73.seqs.biallelic.pm3bp.crossmapped.merged.bed; done
```
(NOTE: The above command is left unaltered for sake of documentation, but does not account for potential multimapping of coordinates. If you would like to repeat this analysis, use the command below instead. In our case, with and without the extra filters did not make a difference.)
```{bash}
for g in [NAMparents...]; do echo "${g}"; gawk -v OFS='\t' '{split($4, a, "_"); chr[$4]=$1; if(a[6]=="St1"){ if(!(a[1]"_"a[2]"_"a[3]"_"a[4]"_"a[5] in st1)){ st1[a[1]"_"a[2]"_"a[3]"_"a[4]"_"a[5]]=$3 } else { multimapped[a[1]"_"a[2]"_"a[3]"_"a[4]"_"a[5]]=a[1]"_"a[2]"_"a[3]"_"a[4]"_"a[5] } } else if(a[6]=="St2"){ if(!(a[1]"_"a[2]"_"a[3]"_"a[4]"_"a[5] in st2)){ st2[a[1]"_"a[2]"_"a[3]"_"a[4]"_"a[5]]=$3 } else { multimapped[a[1]"_"a[2]"_"a[3]"_"a[4]"_"a[5]]=a[1]"_"a[2]"_"a[3]"_"a[4]"_"a[5] } } else if(a[6]=="Sp1"){ if(!(a[1]"_"a[2]"_"a[3]"_"a[4]"_"a[5] in sp1)){ sp1[a[1]"_"a[2]"_"a[3]"_"a[4]"_"a[5]]=$3 } else { multimapped[a[1]"_"a[2]"_"a[3]"_"a[4]"_"a[5]]=a[1]"_"a[2]"_"a[3]"_"a[4]"_"a[5] } } else if(a[6]=="Sp2"){ if(!(a[1]"_"a[2]"_"a[3]"_"a[4]"_"a[5] in sp2)){ sp2[a[1]"_"a[2]"_"a[3]"_"a[4]"_"a[5]]=$3 } else { multimapped[a[1]"_"a[2]"_"a[3]"_"a[4]"_"a[5]]=a[1]"_"a[2]"_"a[3]"_"a[4]"_"a[5] } }} END { for(i in st1){ if((i in st2) && (i in sp1) && (i in sp2)) { if(!(i in multimapped)){ if((chr[i"_St1"]==chr[i"_St2"]) && (chr[i"_St1"]==chr[i"_Sp1"]) && (chr[i"_St1"]==chr[i"_Sp2"])){ if(st1[i]<st2[i]){ print chr[i"_St1"], st1[i], st2[i], sp1[i], sp2[i], i } else if(st1[i]>st2[i]){ print chr[i"_St1"], sp2[i], sp1[i], st2[i], st1[i], i }}}}}}' ${g}_againstB73.seqs.biallelic.pm3bp.crossmapped.bed | sort -k1,1 -k2,2n > ${g}_againstB73.seqs.biallelic.pm3bp.crossmapped.merged_multimapFIXed.bed; done
```


<br/><br/>

If genome carries the insertion allele, take the coordinates of the insertion itself. If genome carries the deletion allele, get +/-3 bp around the deletion site.:

```{bash}
for g in [NAMparents...]; do echo "${g}"; gawk -v OFS='\t' '{ if(($4-$3)==1){ print $1, $2-1, $5, $6"_del" } else if(($4-$3)>1){ print $1, $3, $4-1, $6"_ins" } }' ./${g}/${g}_againstB73.seqs.biallelic.pm3bp.crossmapped.merged.bed | sort -k1,1 -k2,2n > ./${g}/${g}_againstB73.seqs.biallelic.pm3bp.crossmapped.merged.pm3Ins.bed; done
```

<br/><br/>

Check if these sites overlap in the NAM genomes coordinates (and therefore would not be biallelic):

```{bash}
for g in [NAMparents...]; do echo "${g}"; bedtools intersect -a ${g}_againstB73.seqs.biallelic.pm3bp.crossmapped.merged.pm3Ins.bed -b ${g}_againstB73.seqs.biallelic.pm3bp.crossmapped.merged.pm3Ins.bed -wo > ${g}_againstB73.seqs.biallelic.pm3bp.crossmapped.merged.pm3Ins.intersect.bed; done
```

<br/><br/>

From the intersected file, get those that need to be removed in an extra file:

```{bash}
for g in [NAMparents...]; do gawk -v OFS='\t' '{if(!($4==$8)){split($4, a, "_"); split($8, b, "_"); print a[1]"_"a[2]"_"a[3]"_"a[4]"_"a[5] "\n" b[1]"_"b[2]"_"b[3]"_"b[4]"_"b[5] }}' ${g}_againstB73.seqs.biallelic.pm3bp.crossmapped.merged.pm3Ins.intersect.bed | sort | uniq > ${g}_againstB73.seqs.biallelic.pm3bp.crossmapped.merged.pm3Ins.intersect.overlapList.bed; done
```

<br/><br/>

Get list of those, where hybrid shows NAM allele, to check if still >=2 hybrids share this allele:

```{bash}
for g in [NAMparents...]; do gawk -v OFS='\t' '{if(!($4==$8)){print $4 "\n" $8}}' ${g}_againstB73.seqs.biallelic.pm3bp.crossmapped.merged.pm3Ins.intersect.bed | sort | uniq | gawk '{split($1, a, ":|_"); if( ((length(a[2])==1) && (a[7]=="del")) || ((length(a[2])>1) && (a[7]=="ins")) ){ print a[1]":"a[2]"_"a[3]"_"a[4]"_"a[5]"_"a[6] }}' > ${g}_againstB73.seqs.biallelic.pm3bp.crossmapped.merged.pm3Ins.intersect.overlapList.withNAMAllel.bed; done

cat B97_againstB73.seqs.biallelic.pm3bp.crossmapped.merged.pm3Ins.intersect.overlapList.withNAMAllel.bed CML103_againstB73.seqs.biallelic.pm3bp.crossmapped.merged.pm3Ins.intersect.overlapList.withNAMAllel.bed CML247_againstB73.seqs.biallelic.pm3bp.crossmapped.merged.pm3Ins.intersect.overlapList.withNAMAllel.bed CML277_againstB73.seqs.biallelic.pm3bp.crossmapped.merged.pm3Ins.intersect.overlapList.withNAMAllel.bed CML333_againstB73.seqs.biallelic.pm3bp.crossmapped.merged.pm3Ins.intersect.overlapList.withNAMAllel.bed M37W_againstB73.seqs.biallelic.pm3bp.crossmapped.merged.pm3Ins.intersect.overlapList.withNAMAllel.bed Mo18W_againstB73.seqs.biallelic.pm3bp.crossmapped.merged.pm3Ins.intersect.overlapList.withNAMAllel.bed Ms71_againstB73.seqs.biallelic.pm3bp.crossmapped.merged.pm3Ins.intersect.overlapList.withNAMAllel.bed NC358_againstB73.seqs.biallelic.pm3bp.crossmapped.merged.pm3Ins.intersect.overlapList.withNAMAllel.bed Oh43_againstB73.seqs.biallelic.pm3bp.crossmapped.merged.pm3Ins.intersect.overlapList.withNAMAllel.bed Tx303_againstB73.seqs.biallelic.pm3bp.crossmapped.merged.pm3Ins.intersect.overlapList.withNAMAllel.bed CML69_againstB73.seqs.biallelic.pm3bp.crossmapped.merged.pm3Ins.intersect.overlapList.withNAMAllel.bed HP301_againstB73.seqs.biallelic.pm3bp.crossmapped.merged.pm3Ins.intersect.overlapList.withNAMAllel.bed Ki11_againstB73.seqs.biallelic.pm3bp.crossmapped.merged.pm3Ins.intersect.overlapList.withNAMAllel.bed Ki3_againstB73.seqs.biallelic.pm3bp.crossmapped.merged.pm3Ins.intersect.overlapList.withNAMAllel.bed Ky21_againstB73.seqs.biallelic.pm3bp.crossmapped.merged.pm3Ins.intersect.overlapList.withNAMAllel.bed M162W_againstB73.seqs.biallelic.pm3bp.crossmapped.merged.pm3Ins.intersect.overlapList.withNAMAllel.bed CML322_againstB73.seqs.biallelic.pm3bp.crossmapped.merged.pm3Ins.intersect.overlapList.withNAMAllel.bed IL14H_againstB73.seqs.biallelic.pm3bp.crossmapped.merged.pm3Ins.intersect.overlapList.withNAMAllel.bed Oh7b_againstB73.seqs.biallelic.pm3bp.crossmapped.merged.pm3Ins.intersect.overlapList.withNAMAllel.bed P39_againstB73.seqs.biallelic.pm3bp.crossmapped.merged.pm3Ins.intersect.overlapList.withNAMAllel.bed A188_againstB73.seqs.biallelic.pm3bp.crossmapped.merged.pm3Ins.intersect.overlapList.withNAMAllel.bed A619_againstB73.seqs.biallelic.pm3bp.crossmapped.merged.pm3Ins.intersect.overlapList.withNAMAllel.bed Mo17_againstB73.seqs.biallelic.pm3bp.crossmapped.merged.pm3Ins.intersect.overlapList.withNAMAllel.bed W22_againstB73.seqs.biallelic.pm3bp.crossmapped.merged.pm3Ins.intersect.overlapList.withNAMAllel.bed | sort | uniq -c > All.seqs.biallelic.pm3bp.crossmapped.merged.pm3Ins.intersect.overlapList.withNAMAllel.bed


gawk -v OFS='\t' '{split($2, a, "_"); if((a[2]-$1)<2){ print $2 }}' All.seqs.biallelic.pm3bp.crossmapped.merged.pm3Ins.intersect.overlapList.withNAMAllel.bed > All.seqs.biallelic.pm3bp.crossmapped.merged.pm3Ins.intersect.overlapList.withNAMAllel.toRemove.bed
```

<br/><br/>

Create one file for every hybrid with those that have to be removed:

```{bash}
for g in [NAMparents...]; do cat All.seqs.biallelic.pm3bp.crossmapped.merged.pm3Ins.intersect.overlapList.withNAMAllel.toRemove.bed ${g}_againstB73.seqs.biallelic.pm3bp.crossmapped.merged.pm3Ins.intersect.overlapList.bed > ${g}_againstB73.seqs.biallelic.pm3bp.toRemove.bed; done
```

<br/><br/>

Remove sites from the list of every hybrid; create one file for each coordinate system:

```{bash}
for g in [NAMparents...]; do gawk -v OFS='\t' '{if(NR==FNR){rem[$1]=$1; next} split($4, ids, "_"); if(!(ids[1]"_"ids[2]"_"ids[3]"_"ids[4]"_"ids[5] in rem)){ print $1, $2, $3, ids[1]"_"ids[2]"_"ids[3]"_"ids[4]"_"ids[5] }}' ${g}_againstB73.seqs.biallelic.pm3bp.toRemove.bed ${g}_againstB73.seqs.biallelic.pm3bp.crossmapped.merged.pm3Ins.bed | sort -k1,1 -k2,2n > ${g}_againstB73.seqs.biallelic.pm3bp.crossmapped.merged.pm3Ins.filtered.${g}coords.bed; done

for g in [NAMparents...]; do gawk -v OFS='\t' '{split($4, a, ":|_"); if(length(a[1])==1){ print a[4], a[5]-2, a[6]+2, $4} else if(length(a[1])>1){ print a[4], a[5]+1, a[6]-1, $4}}' ${g}_againstB73.seqs.biallelic.pm3bp.crossmapped.merged.pm3Ins.filtered.${g}coords.bed | sort -k1,1 -k2,2n > ${g}_againstB73.seqs.biallelic.pm3bp.crossmapped.merged.pm3Ins.filtered.B73coords.bed; done
```


<br/><br/>

## Step 4: Determine MOA-seq coverage at/around INDELs

NAM alleles do not have correct chromosome names yet and both need to be sorted but not for A188 A619 Mo17 W22:

```{bash}
for g in [NAMparents...]; do gawk -v OFS='\t' -v g=${g} '{print g"-chr" $0}' ${g}_againstB73.seqs.biallelic.pm3bp.crossmapped.merged.pm3Ins.filtered.${g}coords.bed | sortBed -g ref_B73${g}.fasta.size.new.txt > ${g}_againstB73.seqs.biallelic.crossmapped.merged.pm3Ins.filtered.${g}coords.chrm_names.sorted.bed; done

for g in [NAMparents...]; do sortBed -g ref_B73${g}.fasta.size.new.txt  -i ${g}_againstB73.seqs.biallelic.pm3bp.crossmapped.merged.pm3Ins.filtered.B73coords.bed > ${g}_againstB73.seqs.biallelic.crossmapped.merged.pm3Ins.filtered.B73coords.sorted.bed; done
```
NOTE: ref_B73${g}.fasta.size.new.txt is a tab-separated file, containing the contig name in the first column and the contig size in the second.

<br/><br/>

This is followed up by the steps described here: https://github.com/jengelhorn/AS-MOA

Only using the scripts which are adapted to INDELs and provided here.



<br/><br/>

## Step 5: Calculate methylation at/around INDELs

This takes the output of the previous step and performes the steps similarly to the analysis with SNPs (see ../methylationAnalysis).

<br/><br/>

Create BED files from the output of the last step (B73.${g}.${tr}.q255.PF.GT.RN.csv):

```{bash}
for g in [NAMparents...]; do for tr in WW DS; do echo "${g}_${tr}"; gawk -v OFS='\t' '{if(!($7=="het")){ print $1, $2-1, $2+length($3), $11, $1":"$2";"$3";"$4";"$5";"$6 }}' B73.${g}.${tr}.q255.PF.GT.RN.csv | sort -k1,1 -k2,2n > ${g}.${tr}.q255.B73coord.bed; done; done

for g in [NAMparents...]; do for tr in WW DS; do echo "${g}_${tr}"; gawk -v OFS='\t' '{if(!($7=="het")){ split($5, namcoord, ":"); if($6=="0/0"){ print namcoord[1], namcoord[2]-1, namcoord[2]+length($3), $11, $1":"$2";"$3";"$4";"$5";"$6 } else if($6=="1/1"){ print namcoord[1], namcoord[2]-1, namcoord[2]+length($4), $11, $1":"$2";"$3";"$4";"$5";"$6 }}}' B73.${g}.${tr}.q255.PF.GT.RN.csv | sort -k1,1 -k2,2n > ${g}.${tr}.q255.${g}coord.bed; done; done
```
The new colums here are: "Chr", "Start", "Stop", "BF", "REF;ALT;NAMID;GT"
Where "Chr" is the B73 chromosome, "Start" is the B73 start position, "Stop" is the B73 stop position, "BF" is the binding frequency, "REF" stands for the reference/B73 allele's sequence, ALT stands for the alternative/NAM-parent allele's sequence, NAMID are the chromosome and stop coordinate in the coordinate system of the NAM parent, and "GT" stands for genotype, which can be 0/0 (indicating homozygous for the reference allele), or 1/1 (indicating heterozygosity).


<br/><br/>

Create BED files with windows around the sites, for which to calculate the average methylation level. The windows are created by adding 20 bp up- and downstream of the site, for both the insertion and the deletion allele.

```{bash}
for g in [NAMparents...]; do for tr in WW DS; do echo "${g}_${tr}"; gawk -v OFS='\t' '{if($2>20){print $1,$2-20,$3+19,$4,$1":"$3";"$5} else {print $1,0,$3+19,$4,$1":"$3";"$5}}' ${g}.${tr}.q255.B73coord.bed | sortBed -g ref_B73${g}.fasta.size.new.sort.txt > ${g}.${tr}.q255.B73coord.41bp.bed; done; done

for g in [NAMparents...]; do for tr in WW DS; do echo "${g}_${tr}"; gawk -v OFS='\t' '{if($2>20){print $1,$2-20,$3+19,$4,$1":"$3";"$5} else {print $1,0,$3+19,$4,$1":"$3";"$5}}' ${g}.${tr}.q255.${g}coord.bed > ${g}.${tr}.q255.${g}coord.41bp.bed; done; done
```


<br/><br/>

Calculate the average methylation over the given windows (here only the CG context is shown as an example):


```{bash}
for g in [NAMparents...]; do for tr in WW DS; do echo "${g}_${tr}"; intersectBed -a ${g}.${tr}.q255.B73coord.41bp.bed -b B73_CG.diploid.bedgraph -wao | gawk -v OFS='\t' '{if(NR==1){ if($9=="."){CHR=$1;St=$2;Sp=$3;BF=$4;ID=$5;n=1;AddVal=0} else {CHR=$1;St=$2;Sp=$3;BF=$4;ID=$5;n=$10;AddVal=$9*$10} } else { if($5==ID){n=n+$10;AddVal=AddVal+($9*$10)} else {print CHR,St,Sp,BF,ID,AddVal/n; if($9=="."){CHR=$1;St=$2;Sp=$3;BF=$4;ID=$5;n=1;AddVal=0} else {CHR=$1;St=$2;Sp=$3;BF=$4;ID=$5;n=$10;AddVal=$9*$10}}}} END{print CHR,St,Sp,BF,ID,AddVal/n}' > ${g}.${tr}.q255.B73coord.41bp.CG.bed; done; done

for g in [NAMparents...]; do for tr in WW DS; do echo "${g}_${tr}"; intersectBed -a ${g}.${tr}.q255.${g}coord.41bp.bed -b ${g}_CG.diploid.bedgraph -wao | gawk -v OFS='\t' '{if(NR==1){ if($9=="."){CHR=$1;St=$2;Sp=$3;BF=$4;ID=$5;n=1;AddVal=0} else {CHR=$1;St=$2;Sp=$3;BF=$4;ID=$5;n=$10;AddVal=$9*$10} } else { if($5==ID){n=n+$10;AddVal=AddVal+($9*$10)} else {print CHR,St,Sp,BF,ID,AddVal/n; if($9=="."){CHR=$1;St=$2;Sp=$3;BF=$4;ID=$5;n=1;AddVal=0} else {CHR=$1;St=$2;Sp=$3;BF=$4;ID=$5;n=$10;AddVal=$9*$10}}}} END{print CHR,St,Sp,BF,ID,AddVal/n}' > ${g}.${tr}.q255.${g}coord.41bp.CG.bed; done; done
```


<br/><br/>

Join methylation of both the B73 and the NAM allele into one file

The new columns are: "B73_chr", "B73_pos", "B73_allele", "NAM_allele", "NAM_ID", "GT", "AMP", "BF", "B73_meth", "NAM_meth", "meth_diff"

where "B73_chr", "B73_pos" are chromosome and stop position

"B73_allele" and "NAM_allele" are the sequences of each allele

"NAM_ID" is the chromosome and stop position in non-B73 parent coordinates

"GT" is the genotype: 0/0 is homozygous B73 and 1/1 is heterozygous

"AMP" this column is all NA for the moment, a classification of if the site is an AMP or not can be added here later

"BF" is the binding frequency

"B73_meth", "NAM_meth" are the calculated methylation over the given window as it is in the B73 or the non-B73 parent

"meth_diff" is the difference in methylation between the parents over this window, i.e. "B73_meth" - "NAM_meth"

```{bash}
for g in [NAMparents...]; do for tr in WW DS; do echo "${g}_${tr}"; gawk -v OFS='\t' -v g=$g '{if(NR==FNR) {split($5, a, ";"); B73[a[2]]=$6; next} split($5, b, ";|:"); print b[3], b[4], b[5], b[6], b[7]":"b[8], b[9], "NA", $4, B73[b[3]":"b[4]], $6, B73[b[3]":"b[4]]-$6 }' ${g}.${tr}.q255.B73coord.41bp.CG.bed ${g}.${tr}.q255.${g}coord.41bp.CG.bed | sort -k1,1 -k2,2n > ${g}.${tr}.q255.41bp.CG.bed; done; done
```











