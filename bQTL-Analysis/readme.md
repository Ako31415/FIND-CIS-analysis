# bQTL-Analysis
This pipeline describes the steps to identify binding quantitative trait loci (bQTL) from transcription factor (TF) binding data in a population of F1 hybrids via local association mapping. Through the addition of parental DNA methylation data, this pipeline allows for identification of bQTL based on genotype, methylation or both. Linear modelling is used to identify significant associations of TF binding frequqncies (BF) in the F1 hybrids (Binding maternal allele/(Bindinding maternal allele + Binding paternal allele)) with either the genotype or the methylation state at SNP or INDEL positions. The scripts provided here are an example application for a population of maize F1 hybrids with a common mother (B73) and 25 diverse paternal lines for which we generated TF binding data via MOA-seq in well-watered and drought conditions.

## Special requirements
- julia (the scripts provided here were tested on version 1.8.1)
- R 4.4.1 (version number is critical)

## Input data

1. Genotype information

This file provides information on the presence or absence of sequence variants (SNPs or INDELs) in the paternal genomes compared to the maternal reference genome. Only positions with at least 2 of the paternal lines carrying the variant while at the same time having TF binding coverage are listed in these files (e.g. MOA polymorphis or MPs from github.com/jengelhorn/AS-MOA). One file per chromosome is expected in the following format:
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
with the columns denoting the B73(maternal) chromosome, the B73 position, the paternal chromosome and position (not need for the analysis), the genotype (not needed for the analysis, information is drawn from the genotype file, for informaiton only: 1/1 for presence of variant, 0/0 for absence), the BF from the MOA-seq analysis (mapped with allowing unique mapping and mapping equally once to both genomes, n.p. in front of the value indicated that the coverage at this position was not high enough to be called as a peak position or was below our read treshold (see github.com/jengelhorn/AS-MOA)), the total normalised MOA-seq read count (maternal+paternal allele), the difference in avergae CG methylation at Cs in 40 bp sourrounding the variant (Pat-Mat), the difference in average CHG methylation at Cs in 40 bp sourrounding the variant (Pat-Mat) and the difference in avergae CHH methylation at Cs in 40 bp sourrounding the variant (Pat-Mat). For lines where DNA methylation data is not available (A619 in our case), the last 3 columns can be omitted. For descriptions of how methylation is counted within a 40 bp window see the methylationAnlaysis folder of this repository.

## Step 1: Data preparation

Script: 1-moa_windows_RUN.sh
Associated scripts: 1-moa_windows_perpare.jl

The script will integrate genotype data (SNP and/or INDEL presence) with relative MOA BF data (condition dependent output for well-watered and drought data separated) and CG, CHG, and CHH methylation data for each genotype except A619. For A619 no methylation data was available and hence separate scripts were prepared to push A619 data into the output data set.

Variables to be set in 1-moa_windows_RUN.sh: 	
- t = e.g., “WW” or “DS”  (environmental condition e.g., WW (well-watered) or DS (drought))
- v3 = “keep” or “”  (whether MOA data in regions that were not in significant peaks (below peak cutoff) should be used ("NPNRtoValue") or dismissed as NA ("NPNRtoNA"))
- PATH_TO_DATA: a directory containing a subdirectory /$t containing the TF binding and methylation data files described above
- PATH_TO_JULIA: the path to the julia installation

Output: $t-MOA_peak_ratio.csv, $t-CG_ratio.csv, $t-CHG_ratio.csv, $t-CHH_ratio.csv, $t-ReadDepth_ratio.csv (depreciated value for read depth at loci) all in a folder named PATH_TO_DATA/BindingFrequency/NPNRtoValue for keeping below peak values and PATH_TO_DATA/BindingFrequency/NPNRtoValue for NA replacement

## Optional steps for no mehtylation data lines:
Script: 1-zPush_A619_into_ratio.jl ; 1-zPush_A619_into_ReadDepth.jl 

Variables to be set in 1-zPush_A619_into_ratio.jl and 1-zPush_A619_into_ReadDepth.jl
- PATH_TO_A619_DATA_WW/A619_WW.csv: the location of the file containing the A619 BF and read depth data for WW, same needs to be done for drought
- PATH_TO_DATA: the same path as used before to prepare the data with methylation for the other lines
- leave "keep" in if it was used in the previous step, otherwise change to "not_keep"

Output: same files as for before but with A619 data added

These scripts will integrate the genotype data (SNP and/or INDEL values) with the relative MOA peak (condition dependent output for well-watered and drought data separated) for A619 while adding “NA” for the missing methylation data. This script is only needed for lines with missing input data such as methylation. 

## Step 2: Spliting data for parallel processing

Script: 2-split.sh

Input:  $t-MOA_peak_ratio.csv, $t-CG_ratio.csv, $t-CHG_ratio.csv, $t-CHH_ratio.csv, $t-ReadDepth_ratio.csv (depreciated value for read depth at loci) all in a folder named PATH_TO_DATA/BindingFrequency/NPNRtoValue for keeping below peak values and PATH_TO_DATA/BindingFrequency/NPNRtoValue for NA replacement

Variables to be set in 2-split.sh: 	
- v1 = “keep” or “”  (whether MOA data in regions that were not in significant peaks (below peak cutoff) should be used ("NPNRtoValue”) or dismissed as NA (“NPNRtoNA”))
- v2 names of the conditions used, in our cae WW and DS
- v3 "MOA_peak” “CG” “CHG” “CHH” “ReadDepth" (input parameters to split)

Output: Splits the integrated MOA, methylation, and read depth data (depreciated) by chromosome in a separate folder. Since the association is tested locally for each input polymorphism (SNP or INDEL) in the genotype data, chromosomes can split to allow paralleled processing.

## Step 3: Local association mapping
 
These are the scripts for batch submission to run the linear regression model to provide raw pvalues (not adjusted for multiple testing) for associations between the haplotype-specific MOA and/or DNA methylation differences at each polymorphism (SNPs and/or INDELs) provided in the input. Two options are provided to either only test associations with genotype only (only SNPs/INDELs are tested, script_Mety_singleFactor_inLM_Model_3.jl) or also for each methylation type (SNPs/INDELs, CG, CHG, CHH; script_error_removed_methlyation_pvalue_issue_fixed_MMtoLM.jl). We note that the script is only provided for "NPNRtoValue” in the output folder and should be adapted for “NPNRtoNA” if needed.

### For genotype association
Script: 3-parallel_mapping_server_Value.LM_geno.sh

Associated script: 3-LM_geno.jl
### For methylation association
Script: 3-parallel_mapping_server_Value.LM_methylation.sh

Associated script: 3-LM_methylation.jl

Variables to change in each bash script:
- PATH_TO_DATA: same as before
- PATH_TO_JULIA
- k: the conditions, in our case WW and DS
- 1 = {1..n} chromosome ID, here 10 chromosomes for maize

Variables to be changed in the julia scripts:
- PATH_TO_DATA: same as before
  
Input files are the output from the previous scripts:
aom = "PATH_TO_DATA/BindingFrequency/NPNRtoValue/splitted/${k}-MOA_peak_ratio_${i}_file.csv" 

rd = "PATH_TO_DATA/BindingFrequency/NPNRtoValue/splitted/${k}-ReadDepth_ratio_${i}_file.csv" 

oneg = "PATH_TO_DATA/GenotypeData_${k}/genotypes_divided_2FPs_${i}.csv"

Outputs: 
One file per chromosome for the genotype association (PATH_TO_DATA/Results/NPNRtoValue/snponly folder) and methylation assocaition (/Results/NPNRtoValue/LM folder). 
Methylation association files contain the following columns:

MOA_window (the position of the variant analysed),log10ProbCHGsingle (log 10 of the p-value for CHG methylation), log10ProbCGsingle (log 10 of the p-value for CG methylation),log10ProbCHHsingle  (log 10 of the p-value for CHH methylation),log10ProbCGcombi (log 10 of the p-value  for combination CG & genotype, depreciated),log10ProbCHGcombi (log 10 of the p-value  for combination CHG & genotype, depreciated),log10ProbCHHcombi (log 10 of the p-value  for combination CHH & genotype, depreciated),ExpVarianceCHG (Variance explained CHG),ExpVarianceCG  (Variance explained CG),ExpVarianceCHH  (Variance explained CHH),ExpVariacneMix  (Variance explained all, depreciated)

Genotype association files contain the following columns:

MOA_window (the position of the variant analysed),SNP (the position of the variant, repeated for technical reasons),log10ProbGeno (log 10 of the p-value for asspciation with the genotype),numberOFsnps (column for technical reasons, can be ignored),ExpVariance (Variance explained).

These files can then be merged for all chromosomes to perform FDR correction of p-values (e.g. in R).






