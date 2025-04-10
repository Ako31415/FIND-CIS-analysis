#!/bin/bash

# run the julia code for all combinations of files. 

for t in WW DS #WW #
do 

v1="/netscratch/dep_psl/grp_frommer/Thomas/Results/HybMoa_0819_WWvsDS/bQTLs_25/NG_revision/$t/"
for v3 in "keep" #"keep_not" #
do

for j in $v1
do

for l in $v3
do

echo "$j $l $t"


	if [ $v3 = "keep" ]
	then loc="NPNRtoValue"
	else loc="NPNRtoNA"
	fi

	if [ ! -d "/netscratch/dep_psl/grp_frommer/Thomas/Results/HybMoa_0819_WWvsDS/bQTLs_25/NG_revision/BindingFrequency/$loc/" ]
	then mkdir -p "/netscratch/dep_psl/grp_frommer/Thomas/Results/HybMoa_0819_WWvsDS/bQTLs_25/NG_revision/BindingFrequency/$loc/"
	fi


 /netscratch/dep_psl/grp_frommer/Thomas/bin/julia/bin/julia -t 80 /biodata/dep_psl/grp_frommer/MOA_raw/Scripts/MichaelS/25_Lines/NG_revision/1-moa_windows_perpare.jl $j $t $l 

done
done
done
done