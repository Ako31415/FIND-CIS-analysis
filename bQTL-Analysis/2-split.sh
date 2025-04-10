set -uex
# split the geno and moa file by chr and run them in parallel
for v1 in "keep" #"keep_not"
do

	if [ $v1 = "keep" ]
	then loc="NPNRtoValue"
	else loc="NPNRtoNA"
	fi

for v2 in "WW" DS
do

head -n1 /netscratch/dep_psl/grp_frommer/Thomas/Results/HybMoa_0819_WWvsDS/bQTLs_25/NG_revision/BindingFrequency/$loc/$v2-CG_ratio.csv > /netscratch/dep_psl/grp_frommer/Thomas/Results/HybMoa_0819_WWvsDS/bQTLs_25/NG_revision/BindingFrequency/$loc/CG_head.csv
head -n1 /netscratch/dep_psl/grp_frommer/Thomas/Results/HybMoa_0819_WWvsDS/bQTLs_25/NG_revision/BindingFrequency/$loc/$v2-CHG_ratio.csv > /netscratch/dep_psl/grp_frommer/Thomas/Results/HybMoa_0819_WWvsDS/bQTLs_25/NG_revision/BindingFrequency/$loc/CHG_head.csv
head -n1 /netscratch/dep_psl/grp_frommer/Thomas/Results/HybMoa_0819_WWvsDS/bQTLs_25/NG_revision/BindingFrequency/$loc/$v2-CHH_ratio.csv > /netscratch/dep_psl/grp_frommer/Thomas/Results/HybMoa_0819_WWvsDS/bQTLs_25/NG_revision/BindingFrequency/$loc/CHH_head.csv
head -n1 /netscratch/dep_psl/grp_frommer/Thomas/Results/HybMoa_0819_WWvsDS/bQTLs_25/NG_revision/BindingFrequency/$loc/$v2-MOA_peak_ratio.csv > /netscratch/dep_psl/grp_frommer/Thomas/Results/HybMoa_0819_WWvsDS/bQTLs_25/NG_revision/BindingFrequency/$loc/MOA_peak_head.csv
head -n1 /netscratch/dep_psl/grp_frommer/Thomas/Results/HybMoa_0819_WWvsDS/bQTLs_25/NG_revision/BindingFrequency/$loc/$v2-ReadDepth_ratio.csv > /netscratch/dep_psl/grp_frommer/Thomas/Results/HybMoa_0819_WWvsDS/bQTLs_25/NG_revision/BindingFrequency/$loc/ReadDepth_head.csv


v3="MOA_peak CG CHG CHH ReadDepth"

	if [ ! -d "/netscratch/dep_psl/grp_frommer/Thomas/Results/HybMoa_0819_WWvsDS/bQTLs_25/NG_revision/BindingFrequency/$loc/splitted" ]
	then mkdir -p "/netscratch/dep_psl/grp_frommer/Thomas/Results/HybMoa_0819_WWvsDS/bQTLs_25/NG_revision/BindingFrequency/$loc/splitted"
	fi

for i in {1..10}
do

for k in $v2
do

for l in $v3
do

#echo $i $j $k $l
#echo /Data/michael/MOA/Met_final/Normal_res/BindingFrequency/NPNRtoValue${j}-${k}-${l}_ratio.csv

# NPNRtoValue

	if [ $v1 = "keep" ]
	then
	grep chr${i}_ /netscratch/dep_psl/grp_frommer/Thomas/Results/HybMoa_0819_WWvsDS/bQTLs_25/NG_revision/BindingFrequency/$loc/${k}-${l}_ratio.csv > /netscratch/dep_psl/grp_frommer/Thomas/Results/HybMoa_0819_WWvsDS/bQTLs_25/NG_revision/BindingFrequency/$loc/${k}-${l}_ratio_${i}.csv
	cat /netscratch/dep_psl/grp_frommer/Thomas/Results/HybMoa_0819_WWvsDS/bQTLs_25/NG_revision/BindingFrequency/$loc/${l}_head.csv /netscratch/dep_psl/grp_frommer/Thomas/Results/HybMoa_0819_WWvsDS/bQTLs_25/NG_revision/BindingFrequency/$loc/${k}-${l}_ratio_${i}.csv > /netscratch/dep_psl/grp_frommer/Thomas/Results/HybMoa_0819_WWvsDS/bQTLs_25/NG_revision/BindingFrequency/$loc/splitted/${k}-${l}_ratio_${i}_file.csv
	rm /netscratch/dep_psl/grp_frommer/Thomas/Results/HybMoa_0819_WWvsDS/bQTLs_25/NG_revision/BindingFrequency/$loc/${k}-${l}_ratio_${i}.csv
	else
# NPNPRtoNA
	grep chr${i}_ /netscratch/dep_psl/grp_frommer/Thomas/Results/HybMoa_0819_WWvsDS/bQTLs_25/NG_revision/BindingFrequency/NPNRtoNA/${k}-${l}_ratio.csv > /netscratch/dep_psl/grp_frommer/Thomas/Results/HybMoa_0819_WWvsDS/bQTLs_25/NG_revision/BindingFrequency/NPNRtoNA/${k}-${l}_ratio_${i}.csv
	cat /netscratch/dep_psl/grp_frommer/Thomas/Results/HybMoa_0819_WWvsDS/bQTLs_25/NG_revision/BindingFrequency/$loc/${l}_head.csv /netscratch/dep_psl/grp_frommer/Thomas/Results/HybMoa_0819_WWvsDS/bQTLs_25/NG_revision/BindingFrequency/NPNRtoNA/${k}-${l}_ratio_${i}.csv > /netscratch/dep_psl/grp_frommer/Thomas/Results/HybMoa_0819_WWvsDS/bQTLs_25/NG_revision/BindingFrequency/NPNRtoNA/splitted/${k}-${l}_ratio_${i}_file.csv
	rm /netscratch/dep_psl/grp_frommer/Thomas/Results/HybMoa_0819_WWvsDS/bQTLs_25/NG_revision/BindingFrequency/NPNRtoNA/${k}-${l}_ratio_${i}.csv
	fi

done 
done
done

rm /netscratch/dep_psl/grp_frommer/Thomas/Results/HybMoa_0819_WWvsDS/bQTLs_25/NG_revision/BindingFrequency/$loc/CG_head.csv
rm /netscratch/dep_psl/grp_frommer/Thomas/Results/HybMoa_0819_WWvsDS/bQTLs_25/NG_revision/BindingFrequency/$loc/CHG_head.csv
rm /netscratch/dep_psl/grp_frommer/Thomas/Results/HybMoa_0819_WWvsDS/bQTLs_25/NG_revision/BindingFrequency/$loc/CHH_head.csv
rm /netscratch/dep_psl/grp_frommer/Thomas/Results/HybMoa_0819_WWvsDS/bQTLs_25/NG_revision/BindingFrequency/$loc/MOA_peak_head.csv
rm /netscratch/dep_psl/grp_frommer/Thomas/Results/HybMoa_0819_WWvsDS/bQTLs_25/NG_revision/BindingFrequency/$loc/ReadDepth_head.csv

done
done