
# Analysis of methylation at transcription factor binding sites 

Transcription factor (TF) binding sites were determined through MOA-seq and following analysis described here:
https://github.com/jengelhorn/AS-MOA


Binding at TF binding sites (TFBS) was determined quantitatively and allele-specific. This is expressed as the binding frequency (BF), which is defined as (number of B73-allele reads)/(number of B73-allele and NAM reads).

BFs were calculated at SNP positions between B73 and the different NAM parents (different genotypes). Methylation was calculated in a +/- 20 bp window around these SNPs.


The input files for the SNPs were the following:

1) A file containing all biallelic SNPs between the B73 and any NAM parents, which are also 1:1 mappable in at least two of the NAM parents. This file is the one generated as output by the [counts_to_table.sh script](https://github.com/jengelhorn/AS-MOA). 

    Here, we will use the following as name for this file: ${g}_allSNPs.tsv

    (${g} as bash variable for each of the NAM parent/genotype names)

@Julia: Please check that it is these files: 
    
    /netscratch/dep_psl/grp_frommer/Thomas/Results/HybMoa_0819_WWvsDS/custom/counts_STAR_80_EG/${g}/new_25_lines_22_12/B73.${g}.${tr}.q255.PF.GT.RN.csv


2) A file containing a subset of the previous one, which only contains SNPs that were classified as MOA footprint polymorphisms (MPs), i.e. those SNPs located in MOA-peaks and a coverage of >7 RPGC for at least one allele, as well as at least 1 read of both alleles.

    Here, we will use the following as name for these files: ${g}_allMPs.tsv

@Julia: Please check that this is these files: /netscratch/dep_psl/grp_frommer/Thomas/Results/HybMoa_0819_WWvsDS/custom/counts_STAR_80_EG/biases_SNPs/q255_new_12_22/${g}.WW.q255.bino.fdr.CPM7.txt


3) A file that is a further subset of the previous (2nd) one, which only contains those MPs, where binding is significantly biased towards one allele, i.e. allele-specific MOA footprint polymorphisms (AMPs). This is the files created as the final output of the allele-specific MOA-analysis described [here](https://github.com/jengelhorn/AS-MOA). 

    Here, we will use the following as name for these files: ${g}_AMPs.tsv

@Julia: Please check that this is these files: /netscratch/dep_psl/grp_frommer/Thomas/Results/HybMoa_0819_WWvsDS/custom/counts_STAR_80_EG/biases_SNPs/q255_new_12_22/AFPs_new_0724/${g}.AFPs.cleaned.WW.csv


<br/><br/>


For the methylation a bedgraph file is needed for each parental line (B73 and each of the NAM parents) and sequence context.

Here, we will use the following name for these files: ${g}\_methylation\_\${m}.bedgraph

With \${m} representing the sequence context, CG, CHG, CHH.


<br/><br/>


## Prepare bed files with windows around SNPs

1. Make a BED file out of the input files

[NAMparent1 NAMparent2 ...] is a stand-in for a space-separated list of all the NAM parent names.
In our case, alternative chromosomes were marked by containing "het" in their name and are filtered out.
Note that SNPs that are not MPs wrongly get labeled as MP here, but as they are filtered out in the next steps, this does not matter unless you want to use this file as is.

Columns of the outputfile will be:
"Chr", "Start", "Stop", "BF", "REF;ALT;NAMID;GT;AMP"

With REF and ALT being the reference (B73) and alternative allele (the allele that some or all NAM parents share). NAMID contains the coordinates of the SNP in the NAM lines coordinates. GT is the genotype (as defined in the input files) and AMP is the label that says whether a SNP was classified as MP or also as AMP (with the caveat described above).

```{bash}
for g in [NAMparent1 NAMparent2 ...]; do echo ${g}; gawk -v OFS='\t' '{if(NR==FNR){ gsub("\"", ""); AMP[$3"_"$4]=$3"_"$4; next } if(!($7=="het")){ if($1"_"$2 in AMP){print $1, $2-1, $2, $11, $3";"$4";"$5";"$6";AMP"} else {print $1, $2-1, $2, $11, $3";"$4";"$5";"$6";MP"} } }' ${g}_AMPs.tsv ${g}_allSNPs.tsv | sort -k1,1 -k2,2n > ${g}.WW.q255.B73coord.bed; done
```

From these output files the same type of bed file is created with the only difference being that the coordinates are of the NAM parent genome:

```{bash}
for g in [NAMparent1 NAMparent2 ...]; do echo ${g}; gawk -v OFS='\t' '{split($5, a, ";"); split(a[3], b, ":"); print b[1], b[2]-1, b[2], $4, a[1]";"a[2]";"$1":"$3";"a[4]";"a[5]}' ${g}.WW.q255.B73coord.bed | sort -k1,1 -k2,2n > ${g}.WW.q255.${g}coord.bed; done
```


2. Make BED files with +/- 20 bp windows

```{bash}
for g in [NAMparent1 NAMparent2 ...]; do echo ${g}; gawk -v OFS='\t' '{if($2>20){print $1,$2-20,$3+20,$4,$1":"$3";"$5} else {print $1,0,$3+20,$4,$1":"$3";"$5}}' ${g}.WW.q255.B73coord.bed > ${g}.WW.q255.B73coord.41bp.bed; done

for g in [NAMparent1 NAMparent2 ...]; do echo ${g}; gawk -v OFS='\t' '{if($2>20){print $1,$2-20,$3+20,$4,$1":"$3";"$5} else {print $1,0,$3+20,$4,$1":"$3";"$5}}' ${g}.WW.q255.${g}coord.bed > ${g}.WW.q255.${g}coord.41bp.bed; done
```


<br/><br/>


## Calculate mean methylation over the given windows

For this a genome_size_file.txt is required containing the chromosome order; compare [bedtools map documentation](https://bedtools.readthedocs.io/en/latest/content/tools/map.html).

Entries of the bedgraph that are describing the score for multiple bases need to be split into multiple entries for each individual basepair, as the means are calculated wrongly otherwise.


1. Calculate methylation of the B73 allele

```{bash}
for g in [NAMparent1 NAMparent2 ...]; do for m in CG CHG CHH; do echo "${g}_${m}"; gawk -v OFS='\t' '{a=($3 - $2); if(a>1){for(i=a; i>=1; i--){ print $1, ($3 - i), ($3 - i + 1), $4 } } else{ print $0 }}' B73_methylation_${m}.bedgraph | bedtools map -a ${g}.WW.q255.B73coord.41bp.bed -b - -g genome_size_file.txt -c 4 -o mean > ${g}.WW.q255.B73coord.41bp.${m}.bed; done; done
```


2. Calculate methylation of the NAM allele

```{bash}
for g in [NAMparent1 NAMparent2 ...]; do for m in CG CHG CHH; do echo "${g}_${m}"; gawk -v OFS='\t' '{a=($3 - $2); if(a>1){for(i=a; i>=1; i--){ print $1, ($3 - i), ($3 - i + 1), $4 } } else{ print $0 }}' ${g}_methylation_${m}.bedgraph | bedtools map -a ${g}.WW.q255.${g}coord.41bp.bed -b - -g genome_size_file.txt -c 4 -o mean > ${g}.WW.q255.${g}coord.41bp.${m}.bed; done; done
```


3. Finishing touches

Some windows do not contain cytosines in any of the three contexts. For these windows, instead of a value a dot is written in the score column. Here, we replace this with a score of zero, as the window does not have any cytosine methylation.

```{bash}
for m in CG CHG CHH; do for g in [NAMparent1 NAMparent2 ...]; do for coord in B73 ${g}; do echo ${m}"_"${g}"_"${coord}; gawk -v OFS='\t' '{if($6=="."){print $1, $2, $3, $4, $5, "0"} else {print $0}}' ${g}.WW.q255.${coord}coord.41bp.${m}.bed > ${g}.WW.q255.${coord}coord.41bp.${m}.noDots.bed; done; done; done
```

Next, join both the B73 and NAM allele methylation into one file:

```{bash}
for g in [NAMparent1 NAMparent2 ...]; do for m in CG CHG CHH; do echo ${g}"_"${m}; gawk -v OFS='\t' -v g=$g '{if(NR==FNR) {split($5, a, ";"); B73[a[1]]=$6; next} split($5, b, ";|:"); print b[5], b[6], b[3], b[4], b[1]":"b[2], b[7], b[8], $4, B73[b[5]":"b[6]], $6, B73[b[5]":"b[6]]-$6 }' ${g}.WW.q255.B73coord.41bp.${m}.noDots.bed ${g}.WW.q255.${g}coord.41bp.${m}.noDots.bed | sort -k1,1 -k2,2n > ${g}.WW.q255.41bp.${m}.tsv; done; done
```

The columns of the output file are: "B73_chr", "B73_pos", "B73_allele", "NAM_allele", "NAM_ID", "GT", "AMP", "BF", "B73_meth", "NAM_meth", "meth_diff"

Where "B73_chr", "B73_pos" are chromosome and stop position of the SNP;

"B73_allele" and "NAM_allele" are the bases of the SNP for each allele;

"NAM_ID" is the chromosome and stop position in non-B73 parent coordinates;

"GT" is the genotype: 0/0 is homozygous for the B73 allele and 1/1 is heterozygous;

"AMP" is either equal to 'AMP' if the MP is classified as an AMP or equal 'MP' if not classified as AMP;

"BF" is the binding frequency;

"B73_meth", "NAM_meth" are the calculated mean methylation over a 41bp window over the SNP (SNP +/- 20bp) as it is in the B73 or the non-B73 parent;

"meth_diff" is the difference in methylation between the parents over this window, i.e. "B73_meth" - "NAM_meth".



<br/><br/>


## Additional filtering of sites

1. Filter out sites where there were no MOA-seq reads and/or no MOA-peaks (marked by "n.r." or "n.p." instead of a BF).

```{bash}
for g in [NAMparent1 NAMparent2 ...]; do for m in CG CHG CHH; do echo ${g}"_"${m}; gawk -v OFS='\t' -v g=$g '{if(!($8~"n.")){print $0}}' ${g}.WW.q255.41bp.${m}.tsv > ${g}.WW.q255.41bp.${m}.noNP.tsv; done; done
```


2. Filter so that only MPs with a heterozygous SNP are listed (Genotype 1/1).

```{bash}
for g in [NAMparent1 NAMparent2 ...]; do for m in CG CHG CHH; do echo ${g}"_"${m}; gawk -v OFS='\t' -v g=$g '{if($6=="1/1"){print $0}}' ${g}.WW.q255.41bp.${m}.noNP.tsv > ${g}.WW.q255.41bp.${m}.noNP.GT1.tsv; done; done
```


3. Filter out all sites that are not MPs.

```{bash}
for g in [NAMparent1 NAMparent2 ...]; do for m in CG CHG CHH; do echo ${g}"_"${m}; gawk -v OFS='\t' '{if(NR==FNR){ gsub("\"", ""); MP[$3"_"$4]=$3"_"$4; next } if($1"_"$2 in MP){ print $0 }}' ${g}_allMPs.tsv ${g}.WW.q255.41bp.${m}.noNP.GT1.tsv > ${g}.WW.q255.41bp.${m}.noNP.GT1.justMPs.tsv; done; done
```



<br/><br/>


## Plotting of results

1. Plot the number of MPs per MP/AMP/strict-AMP category that are either CG or CHG differentially methylated or both. Strict AMPs were defined as MPs where the BF is <=0.15 or >=0.85

```{bash}
Rscript --vanilla create_propDM_boxplots_CGCHG.R
```


2. Plot the methylation differences for different BF bins.

Mutate the input file to contain binned data:

```{bash}
for g in [NAMparent1 NAMparent2 ...]; do for m in CG CHG CHH; do echo ${g}"_"${m}; Rscript --vanilla compute_methBFdistribution_tables.R ${g} ${m} "0.025"; done; done
```


Combine files for all lines:

```{bash}
for g in [NAMparent1 NAMparent2 ...]; do echo ${g}; gawk -v OFS=',' '{if(NR>1){print $0}}' ${g}.WW.q255.41bp.CG.noNP.GT1.justMPs.dmProportionsPerBFBin.0.025bins.csv >> ALL_lines.WW.q255.41bp.CG.noNP.GT1.justMPs.dmProportionsPerBFBin.0.025bins.csv; done

for g in [NAMparent1 NAMparent2 ...]; do echo ${g}; gawk -v OFS=',' '{if(NR>1){print $0}}' ${g}.WW.q255.41bp.CHG.noNP.GT1.justMPs.dmProportionsPerBFBin.0.025bins.csv >> ALL_lines.WW.q255.41bp.CHG.noNP.GT1.justMPs.dmProportionsPerBFBin.0.025bins.csv; done

for g in [NAMparent1 NAMparent2 ...]; do echo ${g}; gawk -v OFS=',' '{if(NR>1){print $0}}' ${g}.WW.q255.41bp.CHH.noNP.GT1.justMPs.dmProportionsPerBFBin.0.025bins.csv >> ALL_lines.WW.q255.41bp.CHH.noNP.GT1.justMPs.dmProportionsPerBFBin.0.025bins.csv; done
```


Create plot:

```{bash}
for m in CG CHG CHH; do Rscript --vanilla create_methylBFdistribution_plots.R ${m}; done
```




3. Create a plot investigating if consistent AMP classification is dependent on if a site is differentially methylated or not. Only look at sites that vary between differential, and non-differential methylation throughout the population.

Merge files of different lines (only shown for CG context here):

```{bash}
for g in [NAMparent1 NAMparent2 ...]; do gawk -v OFS='\t' '{print $0}' ${g}.WW.q255.41bp.CG.noNP.GT1.justMPs.bed >> ALL_Lines.WW.q255.41bp.CG.noNP.GT1.justMPs.bed; done
```


Plot:

```{bash}
Rscript --vanilla create_bindingBiasConsistency_plots.R
```






## Supplementary analyses


1. Create plots that compare methylation in 41bp and 11bp window with regard to those who show binding to the hypermethylated allele

Repeat the steps described above with slight alteration to create files listing the mean methylation over +/-5 bp windows (11 bp total), instead of +/-20 bp windows (42 bp total).

Files here will be called (only CG context): ${g}.WW.q255.11bp.CG.noNP.GT1.justMPs.tsv


Merge files of both window sizes:

```{bash}
for g in [NAMparent1 NAMparent2 ...]; do echo ${g}; gawk -v OFS='\t' '{if(NR==FNR){ meth11bp[$1"_"$2]=$9"_"$10"_"$11; next } if($1"_"$2 in meth11bp){ split(meth11bp[$1"_"$2], a, "_"); print $0, a[1], a[2], a[3]}}' ${g}.WW.q255.11bp.CG.noNP.GT1.justMPs.tsv ${g}.WW.q255.41bp.CG.noNP.GT1.justMPs.tsv > ${g}.WW.q255.41bpVS11bp.CG.noNP.GT1.justMPs.tsv; done
```


Add columns that label wether a site has binding towards the hypomethylated allele (hypoBind) or the hypermethylated allele (hyperBind) for each of the window sizes:

```{bash}
for g in [NAMparent1 NAMparent2 ...]; do echo ${g}; gawk -v OFS='\t' '{if($7=="AMP"){ if(($8>0.5 && $9>0.7 && $10<0.1) || ($8<0.5 && $10>0.7 && $9<0.1)) { Cat41bp="hyperBind" } else if(($8<0.5 && $9>0.7 && $10<0.1) || ($8>0.5 && $10>0.7 && $9<0.1)) { Cat41bp="hypoBind" } else { Cat41bp="noDM" }; if(($8>0.5 && $12>0.7 && $13<0.1) || ($8<0.5 && $13>0.7 && $12<0.1)) { Cat11bp="hyperBind" } else if(($8<0.5 && $12>0.7 && $13<0.1) || ($8>0.5 && $13>0.7 && $12<0.1)) { Cat11bp="hypoBind" } else { Cat11bp="noDM" }; print $1, $2, $5, $8, $9, $10, $11, $12, $13, $14, Cat41bp, Cat11bp }}' ${g}.WW.q255.41bpVS11bp.CG.noNP.GT1.justMPs.tsv > ${g}.WW.q255.41bpVS11bp.CG.noNP.GT1.justMPs.AMP_methCats.tsv; done
```


Add columns with extra counts:

```{bash}
for g in [NAMparent1 NAMparent2 ...]; do echo ${g}; gawk -v OFS='\t' -v g=${g} 'BEGIN{ AMPCount=0; hyperbindCount=0; toHyper=0; toHypo=0; toNoDm=0 } {AMPCount+=1; if($11=="hyperBind"){ hyperbindCount+=1; if($12=="hyperBind") { toHyper+=1 } else if($12=="hypoBind") { toHypo+=1 } else if($12=="noDM") { toNoDm+=1 }}} END{ print g, AMPCount, hyperbindCount, toHyper, toHypo, toNoDm, hyperbindCount/AMPCount, toHyper/hyperbindCount, toHypo/hyperbindCount, toNoDm/hyperbindCount }' ${g}.WW.q255.41bpVS11bp.CG.noNP.GT1.justFPs.AMP_methCats.tsv >> ALL_Lines.41bpVS11bp.CG.hyperBindAMPChanges.tsv; done
```

New columns are: "NAM_line", "total_AFP_count", "41bp_hyperbindCount", "to_hyper_count", "to_hypo_count", "to_NoDM_count", "41bp_hyperbind_percent", "to_hyper_percent", "to_hypo_percent", "to_NoDM_percent"

"NAM_line" indicate which hybrid the information is about;
"total_AFP_count" how many AMPs the hybrid has in total;
"41bp_hyperbindCount" how many AMPs have binding preferentially to the hypermethylated allele in total (41 bp window);
"to_hyper_count" how many AMPs have binding preferentially to the hypermethylated allele in total (41 bp window  AND 11 bp window);
"to_hypo_count" how many AMPs have binding preferentially to the hypermethylated allele in in the 41 bp window, but binding to the hypomethylated allele in the 11 bp window;
"to_NoDM_count" how many AMPs have binding preferentially to the hypermethylated allele in in the 41 bp window, but are not differentially methylated in the 11 bp window;
"41bp_hyperbind_percent" ratio of "41bp_hyperbindCount" and the total AMP count;
"to_hyper_percent" ratio of "to_hyper_count" out of all hypermethylated-allele-bound sites for the 41 bp window;
"to_hypo_percent" ratio of "to_hypo_count" out of all hypermethylated-allele-bound sites for the 41 bp window;
"to_NoDM_percent" ratio of "to_NoDM_count" out of all hypermethylated-allele-bound sites for the 41 bp window;


