set -uex
# split the geno and moa file by chr and run them in parallel
for v1 in "keep" #"keep_not"
do

	if [ $v1 = "keep" ]
	then loc="NPNRtoValue"
	else loc="NPNRtoNA"
	fi

for v2 in WW DS
do

head -n1 PATH_TO_DATA/BindingFrequency/$loc/$v2-CG_ratio.csv > PATH_TO_DATA/BindingFrequency/$loc/CG_head.csv
head -n1 PATH_TO_DATA/BindingFrequency/$loc/$v2-CHG_ratio.csv > PATH_TO_DATA/BindingFrequency/$loc/CHG_head.csv
head -n1 PATH_TO_DATA/BindingFrequency/$loc/$v2-CHH_ratio.csv > PATH_TO_DATA/BindingFrequency/$loc/CHH_head.csv
head -n1 PATH_TO_DATA/BindingFrequency/$loc/$v2-MOA_peak_ratio.csv > PATH_TO_DATA/BindingFrequency/$loc/MOA_peak_head.csv
head -n1 PATH_TO_DATA/BindingFrequency/$loc/$v2-ReadDepth_ratio.csv > PATH_TO_DATA/BindingFrequency/$loc/ReadDepth_head.csv


v3="MOA_peak CG CHG CHH ReadDepth"

	if [ ! -d "PATH_TO_DATA/BindingFrequency/$loc/splitted" ]
	then mkdir -p "PATH_TO_DATA/BindingFrequency/$loc/splitted"
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
	grep chr${i}_ PATH_TO_DATA/BindingFrequency/$loc/${k}-${l}_ratio.csv > PATH_TO_DATA/BindingFrequency/$loc/${k}-${l}_ratio_${i}.csv
	cat PATH_TO_DATA/BindingFrequency/$loc/${l}_head.csv PATH_TO_DATA/BindingFrequency/$loc/${k}-${l}_ratio_${i}.csv > PATH_TO_DATA/BindingFrequency/$loc/splitted/${k}-${l}_ratio_${i}_file.csv
	rm PATH_TO_DATA/BindingFrequency/$loc/${k}-${l}_ratio_${i}.csv
	else
# NPNPRtoNA
	grep chr${i}_ PATH_TO_DATA/BindingFrequency/NPNRtoNA/${k}-${l}_ratio.csv > PATH_TO_DATA/BindingFrequency/NPNRtoNA/${k}-${l}_ratio_${i}.csv
	cat PATH_TO_DATA/BindingFrequency/$loc/${l}_head.csv PATH_TO_DATA/BindingFrequency/NPNRtoNA/${k}-${l}_ratio_${i}.csv > PATH_TO_DATA/BindingFrequency/NPNRtoNA/splitted/${k}-${l}_ratio_${i}_file.csv
	rm PATH_TO_DATA/BindingFrequency/NPNRtoNA/${k}-${l}_ratio_${i}.csv
	fi

done 
done
done

rm PATH_TO_DATA/BindingFrequency/$loc/CG_head.csv
rm PATH_TO_DATA/BindingFrequency/$loc/CHG_head.csv
rm PATH_TO_DATA/BindingFrequency/$loc/CHH_head.csv
rm PATH_TO_DATA/BindingFrequency/$loc/MOA_peak_head.csv
rm PATH_TO_DATA/BindingFrequency/$loc/ReadDepth_head.csv

done
done
