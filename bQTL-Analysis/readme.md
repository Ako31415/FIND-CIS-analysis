# bQTL-Analysis
This pipeline describes the steps to identify binging quantitative trait loci (bQTL) from transcription factor (TF) binding data in a population of F1 hybrids via local association mapping. Through the addition of parental DNA methylation data, this pipeline allows for identification of bQTL based on genotype, methylation or both. Linear modelling is used to identify significant associations of TF binding frequqncies in the F1 hybrids(Binding maternal allele/(Bindinding maternal allele + Binding paternal allele)) with either the genotype or the methylation state at SNP or INDEL positions. The scripts provided here are an example application for a population of maize F1 hybrids with a common mother (B73) and 25 diverse paternal lines for which we generated TF binding Data via MOA-seq.

## Input data
1. Genotype information
   Genotype information should be provided in the following format:
`#CHROM,POS,A188,A619,B97,CML103,CML247,CML277,CML322,CML333,CML69,HP301,IL14H,Ki11,Ki3,Ky21,M162W,M37W,Mo17,Mo18W,Ms71,NC358,Oh43,Oh7b,P39,Tx303,W22
B73-chr1,28762,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
B73-chr1,39898,0,-1,0,0,-1,0,-1,-1,-1,0,-1,0,-1,-1,0,0,0,-1,0,0,-1,0,0,0,-1
B73-chr1,124936,-1,0,-1,-1,-1,-1,-1,0,-1,-1,0,-1,-1,-1,-1,-1,-1,-1,-1,-1,0,-1,-1,-1,-1
B73-chr1,124937,-1,-1,-1,-1,0,-1,0,-1,-1,-1,-1,-1,0,-1,-1,-1,-1,0,-1,-1,-1,-1,-1,-1,-1
B73-chr1,124958,-1,-1,0,0,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,0,-1,-1,-1`



3. TF binding and methylation data

## Step 1: Data preparation

Script: 1-moa_windows_RUN.sh
Associated scripts: 1-moa_windows_perpare.jl

The script will integrate genotype data (SNP and/or INDEL presence) with relative MOA peak (condition dependent output for well-watered and drought data separated) and CG, CHG, and CHH methylation data for each genotype except A619.  For A619 no methylation data was available and hence separate scripts were prepared to push A619 data into the output data set.



Input: 

Variables: 	$ENV = e.g., “WW” or “DS”  (environmental condition e.g., WW (well-watered) or DS (drought))
$miss = “keep” or “”  (whether MOA data in regions that were not in significant peaks (below peak cutoff) should be used ("NPNRtoValue) or dismissed as NA (NPNRtoNA))

Output: $ENV-MOA_peak_ratio.csv, $ENV-CG_ratio.csv, $ENV-CHG_ratio.csv, $ENV-CHH_ratio.csv, $ENV-ReadDepth_ratio.csv (depreciated value for read depth at loci)



