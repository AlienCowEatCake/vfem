#!/bin/bash

for real_air in `LANG=C seq 0 1 7`
do
	for imag_air in `LANG=C seq 0 1 7`
	do
		echo -e "${real_air} ${imag_air}\n" > chi_air.txt
		for real_water in `LANG=C seq 0 1 7`
		do
			for imag_water in `LANG=C seq 0 1 7`
			do
				echo -e "${real_water} ${imag_water}\n" > chi_water.txt
				echo "Testing in # Air:(${real_air} , ${imag_air}), Water:(${real_water} , ${imag_water}) ..."
				./vfem | tee "A_${real_air}_${imag_air}_W_${real_water}_${imag_water}.txt"
				pltname=`ls source_pml*.plt | sort | tail -1`
				pltname_new=`echo ${pltname} | sed "s/\./\.A_${real_air}_${imag_air}_W_${real_water}_${imag_water}\./g"`
				mv ${pltname} ${pltname_new}
				datname=`ls source_pml*.dat | sort | tail -1`
				datname_new=`echo ${datname} | sed "s/\./\.A_${real_air}_${imag_air}_W_${real_water}_${imag_water}\./g"`
				mv ${datname} ${datname_new}
			done
		done
	done
done
