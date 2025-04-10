#!/bin/bash

#set -uex

#v1="uni 2AFPs 3AFPs"

for k in WW DS
do

#source make sure R 4.4.1 is installed and the version that is used


# if [ ! -d PATH_TO_DATA/Results/NPNRtoValue/LM/$k ]
# then mkdir -p PATH_TO_DATA/Results/NPNRtoValue/LM/$k
# fi

if [ ! -d PATH_TO_DATA/Results/NPNRtoValue/snponly/$k ]
then mkdir -p PATH_TO_DATA/Results/NPNRtoValue/snponly/$k
fi

for i in {1..10}
do
	# for sp in 10 #{1..20}
	# do
aom="PATH_TO_DATA/BindingFrequency/NPNRtoValue/splitted/${k}-MOA_peak_ratio_${i}_file.csv" 
rd="PATH_TO_DATA/BindingFrequency/NPNRtoValue/splitted/${k}-ReadDepth_ratio_${i}_file.csv" 
oneg="PATH_TO_DATA/GenotypeData_${k}/genotypes_divided_2FPs_${i}.csv"

    
				PATH_TO_JULIA/julia -t 1 \
				3-LM.jl \
				$aom $i $k $rd $oneg & 
done
done
