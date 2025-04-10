# bQTL-Analysis
This pipeline describes the steps to identify binging quantitative trait loci (bQTL) from transcription factor (TF) binding data in a population of F1 hybrids via local association mapping. Through the addition of parental DNA methylation data, this pipeline allows for identification of bQTL based on genotype, methylation or both. Linear modelling is used to identify significant associations of TF binding frequqncies in the F1 hybrids(Binding maternal allele/(Bindinding maternal allele + Binding paternal allele)) with either the genotype or the methylation state at SNP or INDEL positions. The scripts provided here are an example application for a population of maize F1 hybrids with a common mother (B73) and 25 diverse paternal lines for which we generated TF binding Data via MOA-seq in well-watered and drought conditions.

## Input data

1. Genotype information

This file provides information on the presence ot absence of sequence variants (SNPs or INDELs) in the paternal genomes compared to the maternal reference genome. Only positions with at least 2 of the paternal lines carrying the variant while at the same time having TF binding coverage are listed in these files (e.g. MOA polymorphis or MPs from github.com/jengelhorn/AS-MOA). One file per chromosome is expected in the following format:
 ``` 
#CHROM,POS,A188,A619,B97,CML103,CML247,CML277,CML322,CML333,CML69,HP301,IL14H,Ki11,Ki3,Ky21,M162W,M37W,Mo17,Mo18W,Ms71,NC358,Oh43,Oh7b,P39,Tx303,W22
B73-chr1,28762,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
B73-chr1,39898,0,-1,0,0,-1,0,-1,-1,-1,0,-1,0,-1,-1,0,0,0,-1,0,0,-1,0,0,0,-1
B73-chr1,124936,-1,0,-1,-1,-1,-1,-1,0,-1,-1,0,-1,-1,-1,-1,-1,-1,-1,-1,-1,0,-1,-1,-1,-1
B73-chr1,124937,-1,-1,-1,-1,0,-1,0,-1,-1,-1,-1,-1,0,-1,-1,-1,-1,0,-1,-1,-1,-1,-1,-1,-1
B73-chr1,124958,-1,-1,0,0,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,0,-1,-1,-1`
 ``` 
Each column represents one paternal line with 0 denoting presence of a variant while -1 denotes the same genotype as B73 for the respective position.

2. TF binding and methylation data

For each F1, a file containing TF binding data and information of parental DNA methylation is needed in the following format:
``` 
B73-chr1     28762    Oh43-chr1:34921       1/1     0.266666        4.52457 -0.925926       -0.694445       -0.008547
B73-chr1     39898    Oh43-chr1:46057       0/0     n.p.0.529412    2.56392 0       0       -0.00512821
B73-chr1    124936    Oh43-chr1:124474      1/1     0.528572        10.5573 0       0       0
```
with the columns denoting the B73(maternal) chromosome, the B73 position, the paternal chromosome and position (not need for the analysis), the genotype (not needed for the analysis, information is drawn from the genotype file, for informaiton only: 1/1 for presence of variant, 0/0 for absence), the BF from the MOA-seq analysis (mapped with allowing unique mapping and mapping equally once to both genomes), the total normalised MOA-seq read count (maternal+paternal allele), the difference in avergae CG methylation at Cs in 40 bp sourrounding the variant (Pat-Mat), the difference in avergae CHG methylation at Cs in 40 bp sourrounding the variant (Pat-Mat) and the difference in avergae CHH methylation at Cs in 40 bp sourrounding the variant (Pat-Mat). For lines where DNA methylation data is not available (A619 in our case), the last 3 columns can be omitted. 

## Step 1: Data preparation

Script: 1-moa_windows_RUN.sh
Associated scripts: 1-moa_windows_perpare.jl

The script will integrate genotype data (SNP and/or INDEL presence) with relative MOA BF data (condition dependent output for well-watered and drought data separated) and CG, CHG, and CHH methylation data for each genotype except A619. For A619 no methylation data was available and hence separate scripts were prepared to push A619 data into the output data set.

Input: 

Variables: 	$ENV = e.g., “WW” or “DS”  (environmental condition e.g., WW (well-watered) or DS (drought))
$miss = “keep” or “”  (whether MOA data in regions that were not in significant peaks (below peak cutoff) should be used ("NPNRtoValue) or dismissed as NA (NPNRtoNA))

Output: $ENV-MOA_peak_ratio.csv, $ENV-CG_ratio.csv, $ENV-CHG_ratio.csv, $ENV-CHH_ratio.csv, $ENV-ReadDepth_ratio.csv (depreciated value for read depth at loci)

## Optional steps for no mehtylation data lines:
Script: 1-zPush_A619_into_ratio.jl ; 1-zPush_A619_into_ReadDepth.jl 

Output: WW-MOA_peak_ratio.csv, WW -CG_ratio.csv, WW -CHG_ratio.csv, WW -CHH_ratio.csv, WW -ReadDepth_ratio.csv (depreciated value for read depth at loci), DS-MOA_peak_ratio.csv, DS -CG_ratio.csv, DS -CHG_ratio.csv, DS -CHH_ratio.csv, DS -ReadDepth_ratio.csv (depreciated value for read depth at loci)

The scripts will integrate the genotype data (SNP and/or INDEL values) with the relative MOA peak (condition dependent output for well-watered and drought data separated) for A619 while adding “NA” for the missing methylation data. This script is only needed for lines with missing input data such as methylation. 

## Step 2: Spliting data for parallel processing

