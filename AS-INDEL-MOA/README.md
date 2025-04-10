# Analysis pipeline for allele-specific MOA-seq analysis using INDELs instead of SNPs

This part of the pipeline shows how a list of biallelic INDELs is created pairwise genome alignments.

After mapping of MOA-seq reads to the genomes an analysis of allele-specific transcription factor (TF) binding can be performed similarly as for SNPs (described here: https://github.com/jengelhorn/AS-MOA ).

Additionally an analysis of allele-specific methylation can be performed (see ../methylationAnalysis) and sites can also be used for bQTL mapping (see ../bQTL-Analysis).



## Step 1: pairwise genome alignments

Genomes of each NAM parent are aligned to one reference genome (here: B73) using Anchorwave v1.2.2.

1. Gene sequences of B73 are extracted using anchorwave gff2seq

```{bash}
anchorwave gff2seq -r Zm-B73-REFERENCE-NAM-5.0.B73.fa -i Zea_mays.Zm-B73-REFERENCE-NAM-5.0.57.B73-chr.gff3 -o B73_anchors_cds.fa
```

2. Gene sequences are aligned to the B73/reference sequence, and each of the other NAM parents genomes to localise anchors

```{bash}
minimap2 -x splice -t 50 -k 12 -a -p 0.4 -N 20 Zm-B73-REFERENCE-NAM-5.0.B73.fa B73_anchors_cds.fa > B73cds_against_B73.sam

for g in [NAMparents...]; do echo "${g}"; mkdir ${g}; minimap2 -x splice -t 10 -k 12 -a -p 0.4 -N 20 ${g}.pseudomolecules-v1-sm.fasta B73_anchors_cds.fa > B73cds_against_${g}.sam; done
```

3. Pairwise alignments are created using anchorwave proali

```{bash}
for g in [NAMparents...]; do echo "${g}"; anchorwave proali -i Zea_mays.Zm-B73-REFERENCE-NAM-5.0.57.B73-chr.gff3 -as B73_anchors_cds.fa -r Zm-B73-REFERENCE-NAM-5.0.B73.fa -a B73cds_against_${g}.sam -ar B73cds_against_B73.sam -s ${g}.pseudomolecules-v2.1-sm.fasta -n anchors -R 1 -Q 1 -o ${g}_againstB73.maf -t 14; done
```


## Step 2: Variant calling

Variants are extracted from the MAF (multiple alignment file) output files using wgatools wgatools 0.1.0.

```{bash}
for g in [NAMparents...]; do wgatools call ${g}_againstB73.maf -n ${g} -s -l 1 -t 80 -v > ${g}_againstB73.vcf; done
```

Chain files are created to enable lift over of coordinates between genomes:

```{bash}
for g in [NAMparents...]; do wgatools maf2chain ./${g}/${g}_againstB73.maf -t 55 -v > ./${g}/${g}_againstB73.chain; done
```


## Step 3: Determining biallelic INDELs

Create a BED file with only INDELs from the VCF files created in the previous step. Include the base before and after the indel, so as to have two bases on the deletion allele to serve as point of comparison.

```{bash}
for g in [NAMparents...]; do gawk -v OFS='\t' '{if(!($0 ~ "^#")){if((length($4)>1 && length($4)<51) || (length($5)>1 && length($5)<51)){print $1, $2-1, $2+length($4), $1":"$2":"$4":"$5}}}' ${g}_againstB73.vcf > ${g}_againstB73.INDELs.bed; done
```

Biallelic sites were determined as sites where every inbred line with the insertion allele had the same insertion sequence and position. The other allele remaining then was a share deletion allele.

Create one file with the sequences of INDELs as 4th and 5th column:

```{bash}
for g in [NAMparents...]; do gawk -v OFS='\t' '{split($4, a, ":"); print $1, $2, $3, a[3]":"a[4]}' ${g}_againstB73.INDELs.bed >> All_againstB73.INDELs.seqs.tsv; done
```


Use bash sort and uniq -c to find out how often an INDEL occurs between B73 and each of the different NAM inbred lines. If an INDEL varies in sequence at a position there will be two separate entries.

```{bash}
sort All_againstB73.INDELs.seqs.tsv | uniq -c | gawk -v OFS='\t' '{print $1, $2, $3, $4, $5}' > All_againstB73.INDELs.seqs.uniq.tsv
```

Remove the count column in order to run bedtools intersect:

```{bash}
gawk -v OFS='\t' '{print $2, $3, $4, $5}' All_againstB73.INDELs.seqs.uniq.tsv | sort -k1,1 -k2,2n > All_againstB73.INDELs.seqs.uniq.noCounts.bed
```


Run bedtools intersect in order to identify INDELs overlapping in their B73 coordinates to exclude them:

```{bash}
bedtools intersect -a All_againstB73.INDELs.seqs.uniq.noCounts.bed -b All_againstB73.INDELs.seqs.uniq.noCounts.bed -wa -wb -loj > All_againstB73.INDELs.seqs.uniq.noCounts.intersect.tsv
```


Only keep entries that occur only once as overlapping with themselves and remove the rest: 

```{bash}
gawk -v OFS='\t' 'BEGIN{vars["-"]="-"; dups["-"]="-"}{if($1":"$2":"$3":"$4 in vars){ dups[$1":"$2":"$3":"$4]=$1":"$2":"$3":"$4 } else if($1":"$2":"$3":"$4==$5":"$6":"$7":"$8 || $5=="."){ vars[$1":"$2":"$3":"$4]=$1":"$2":"$3":"$4 } else { vars[$1":"$2":"$3":"$4]=$1":"$2":"$3":"$4 }} END{ for(i in vars){ if(!(i in dups)){ split(vars[i], a, ":"); print a[1], a[2], a[3], a[4]":"a[5] }}}' All_againstB73.INDELs.seqs.uniq.noCounts.intersect.tsv | sort -k1,1 -k2,2n > All_againstB73.INDELs.seqs.biallelic.bed
```

Take the counts of how often the INDELs occur in the different hybrids and add them to the biallic file in order to filter for a minimum of 2 hybrids the variants occur in:

```{bash}
gawk -v OFS='\t' '{if(NR==FNR){ counts[$2"_"$3"_"$4"_"$5]=$1; next} if($1"_"$2"_"$3"_"$4 in counts){ print $0, counts[$1"_"$2"_"$3"_"$4] }}' All_againstB73.INDELs.seqs.uniq.bed All_againstB73.INDELs.seqs.biallelic.bed | sort -k1,1 -k2,2n > All_againstB73.INDELs.seqs.biallelic.counts.tsv
```


Filter out positions that only occur in 1 hybrid:

```{bash}
gawk -v OFS='\t' '{if($5>1){ print $1, $2, $3, $4"_"$5 }}' All_againstB73.INDELs.seqs.biallelic.counts.tsv > All_againstB73.INDELs.seqs.biallelic.counts.min2hyb.bed
```

Create BED files with three positions before and after the INDEL. This is so that there are some bases to compare coverage for the deletion allele. Start and Stop coordinates are split into 4 separate rows to make liftover easier:

```{bash}
gawk -v OFS='\t' '{print $1, $2-2, $2-1, $4"_"$1"_"$2"_"$3"_St1\n" $1, $2, $2+1, $4"_"$1"_"$2"_"$3"_St2\n" $1, $3-1, $3, $4"_"$1"_"$2"_"$3"_Sp1\n" $1, $3+1, $3+2, $4"_"$1"_"$2"_"$3"_Sp2"}' All_againstB73.INDELs.seqs.biallelic.counts.min2hyb.bed > All_againstB73.INDELs.seqs.biallelic.counts.min2hyb.pm3bp.bed
```


Lift over files with start and stop positions of the selected INDELs to each NAM parent genome using CrossMap (version 0.7.0). The chain files created in step 2 are used.

```{bash}
for g in [NAMparents...]; do CrossMap bed --chromid a --unmap-file ${g}_againstB73.seqs.biallelic.pm3bp.unmapped.bed ${g}_againstB73.chain All_againstB73.INDELs.seqs.biallelic.counts.min2hyb.pm3bp.bed ${g}_againstB73.seqs.biallelic.pm3bp.crossmapped.bed; done
```

Merge start and stop positions accordingly after lift over, if all 4 coordinates are present:

```{bash}
for g in [NAMparents...]; do echo "${g}"; gawk -v OFS='\t' '{split($4, a, "_"); chr[a[1]"_"a[2]"_"a[3]"_"a[4]"_"a[5]]=$1; if(a[6]=="St1"){ st1[a[1]"_"a[2]"_"a[3]"_"a[4]"_"a[5]]=$3 } else if(a[6]=="St2"){ st2[a[1]"_"a[2]"_"a[3]"_"a[4]"_"a[5]]=$3 } else if(a[6]=="Sp1"){ sp1[a[1]"_"a[2]"_"a[3]"_"a[4]"_"a[5]]=$3 } else if(a[6]=="Sp2"){ sp2[a[1]"_"a[2]"_"a[3]"_"a[4]"_"a[5]]=$3 }} END { for(i in st1){ if((i in st2) && (i in sp1) && (i in sp2)) { if(st1[i]<st2[i]){ print chr[i], st1[i], st2[i], sp1[i], sp2[i], i } else if(st1[i]>st2[i]){ print chr[i], sp2[i], sp1[i], st2[i], st1[i], i }}}}' ${g}_againstB73.seqs.biallelic.pm3bp.crossmapped.bed | sort -k1,1 -k2,2n > ${g}_againstB73.seqs.biallelic.pm3bp.crossmapped.merged.bed; done
```

(If there is an insertion, take the coordinates 

```{bash}
for g in B97 CML103 CML247 CML277 CML333 M37W Mo18W Ms71 NC358 Oh43 Tx303 CML69 HP301 Ki11 Ki3 Ky21 M162W CML322 IL14H Oh7b P39 A188 A619 Mo17 W22; do echo "${g}"; gawk -v OFS='\t' '{ if(($4-$3)==1){ print $1, $2-1, $5, $6"_del" } else if(($4-$3)>1){ print $1, $3, $4-1, $6"_ins" } }' ./${g}/${g}_againstB73.seqs.biallelic.pm3bp.crossmapped.merged.bed | sort -k1,1 -k2,2n > ./${g}/${g}_againstB73.seqs.biallelic.pm3bp.crossmapped.merged.pm3Ins.bed; done
```





















