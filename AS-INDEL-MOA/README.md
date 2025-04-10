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

























