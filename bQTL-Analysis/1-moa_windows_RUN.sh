#!/bin/bash

# run the julia code for all combinations of files. 

for t in WW DS #add names of folders containing files for different conditions here 
do 

v1="PATH_TO_DATA/$t/"
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

	if [ ! -d "PATH_TO_DATA/BindingFrequency/$loc/" ]
	then mkdir -p "PATH_TO_DATA/BindingFrequency/$loc/"
	fi


 PATH_TO_JULIA/bin/julia -t 80 1-moa_windows_perpare.jl $j $t $l 

done
done
done
done
