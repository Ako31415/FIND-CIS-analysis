#!/bin/bash

#set -uex

#v1="uni 2AFPs 3AFPs"

for k in WW DS
do

#source /opt/share/software/scs/appStore/modules/init/profile.sh
source /netscratch/common/Saurabh/envMod/init/bash
module load mambaforge/auto-managed/v23.1.0-3
conda activate /opt/share/software/scs/appStore/selfStorage/thartwig/bookworm/mambaforge/v23.1.0-3/envs/R4.4.1



# if [ ! -d /netscratch/dep_psl/grp_frommer/Thomas/Results/HybMoa_0819_WWvsDS/bQTLs_25/NG_revision/Results/NPNRtoValue/LM/$k ]
# then mkdir -p /netscratch/dep_psl/grp_frommer/Thomas/Results/HybMoa_0819_WWvsDS/bQTLs_25/NG_revision/Results/NPNRtoValue/LM/$k
# fi

if [ ! -d /netscratch/dep_psl/grp_frommer/Thomas/Results/HybMoa_0819_WWvsDS/bQTLs_25/NG_revision/Results/NPNRtoValue/snponly/$k ]
then mkdir -p /netscratch/dep_psl/grp_frommer/Thomas/Results/HybMoa_0819_WWvsDS/bQTLs_25/NG_revision/Results/NPNRtoValue/snponly/$k
fi

for i in {1..10}
do
	# for sp in 10 #{1..20}
	# do
aom="/netscratch/dep_psl/grp_frommer/Thomas/Results/HybMoa_0819_WWvsDS/bQTLs_25/NG_revision/BindingFrequency/NPNRtoValue/splitted/${k}-MOA_peak_ratio_${i}_file.csv" 
rd="/netscratch/dep_psl/grp_frommer/Thomas/Results/HybMoa_0819_WWvsDS/bQTLs_25/NG_revision/BindingFrequency/NPNRtoValue/splitted/${k}-ReadDepth_ratio_${i}_file.csv" 
oneg="/netscratch/dep_psl/grp_frommer/Thomas/Results/HybMoa_0819_WWvsDS/bQTLs_25/NG_revision/GenotypeData_${k}/genotypes_divided_2FPs_${i}.csv"

			LD_PRELOAD=/opt/share/software/scs/appStore/selfStorage/thartwig/bookworm/mambaforge/v23.1.0-3/envs/R4.4.1/x86_64-conda-linux-gnu/lib/libstdc++.so.6 \
				/netscratch/dep_psl/grp_frommer/Thomas/bin/julia/bin/julia -t 1 \
				/biodata/dep_psl/grp_frommer/MOA_raw/Scripts/MichaelS/25_Lines/NG_revision/script_error_removed_methlyation_pvalue_issue_fixed_MMtoLM.jl \
				$aom $i $k $rd $oneg & 
done
done
#NG_revision/script_Mety_singleFactor_inLM_Model_3.jl
		