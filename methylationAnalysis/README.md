
# Analysis of methylation at transcription factor binding sites 

Transcription factor (TF) binding sites were determined through MOA-seq and following analysis described here:
https://github.com/jengelhorn/AS-MOA


Binding at TF binding sites (TFBS) was determined quantitatively and allele-specific. This is expressed as the binding frequency (BF), which is defined as (number of B73-allele reads)/(number of B73-allele and NAM reads).

BFs were calculated at SNP positions between B73 and the different NAM parents (different genotypes). Methylation was calculated in a +/- 20 bp window around these SNPs.


The input files for the SNPs were the following:

1) A file containing all biallelic SNPs between the B73 and any NAM parents, which are also 1:1 mappable in at least two of the NAM parents. This file is the one generated as output by the [counts_to_table.sh script](https://github.com/jengelhorn/AS-MOA). 

    Here, we will use the following as name for this file: ${g}_allSNPs.tsv
(${g} as bash variable for each of the NAM parent/genotype names)

@Julia: Please check that it is these files: /netscratch/dep_psl/grp_frommer/Thomas/Results/HybMoa_0819_WWvsDS/custom/counts_STAR_80_EG/${g}/new_25_lines_22_12/B73.${g}.${tr}.q255.PF.GT.RN.csv


2) A file containing a subset of the previous one, which only contains SNPs that were classified as MOA footprint polymorphisms (MPs), i.e. those SNPs located in MOA-peaks and a coverage of >7 RPGC for at least one allele, as well as at least 1 read of both alleles.

    Here, we will use the following as name for these files: ${g}_allMPs.tsv

@Julia: Please check that this is these files: /netscratch/dep_psl/grp_frommer/Thomas/Results/HybMoa_0819_WWvsDS/custom/counts_STAR_80_EG/biases_SNPs/q255_new_12_22/${g}.WW.q255.bino.fdr.CPM7.txt


3) A file that is a further subset of the previous (2nd) one, which only contains those MPs, where binding is significantly biased towards one allele, i.e. allele-specific MOA footprint polymorphisms (AMPs). This is the files created as the final output of the allele-specific MOA-analysis described [here](https://github.com/jengelhorn/AS-MOA). 

    Here, we will use the following as name for these files: ${g}_AMPs.tsv

@Julia: Please check that this is these files: /netscratch/dep_psl/grp_frommer/Thomas/Results/HybMoa_0819_WWvsDS/custom/counts_STAR_80_EG/biases_SNPs/q255_new_12_22/AFPs_new_0724/${g}.AFPs.cleaned.WW.csv


<br/><br/>


For the methylation a bedgraph file is needed for each parental line (B73 and each of the NAM parents).
Here, we will use the following name for these files: ${g}_methylation.bedgraph


<br/><br/>


## Prepare bed files with windows around SNPs

```{bash}
for g in [NAMparent1 NAMparent2 ...]; do echo ${g}; gawk -v OFS='\t' '{if(NR==FNR){ gsub("\"", ""); AFP[$3"_"$4]=$3"_"$4; next } if(!($7=="het")){ if($1"_"$2 in AFP){print $1, $2-1, $2, $11, $3";"$4";"$5";"$6";AFP"} else {print $1, $2-1, $2, $11, $3";"$4";"$5";"$6";FP"} } }' ${g}_AMPs.tsv ${g}_allSNPs.tsv | sort -k1,1 -k2,2n > ${g}.WW.q255.B73coord.bed; done


```



