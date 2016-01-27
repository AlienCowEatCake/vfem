#!/bin/bash

./vfem_classic | tee "classic.txt"
for real in `LANG=C seq 3.5 0.2 4.5`
do
	for imag in `LANG=C seq 2.5 0.2 3.5`
	do
		echo "Testing in (${real} , ${imag}) ..."
		echo -e "${real} ${imag}\n" > chi.txt
#		./vfem | tee "${real}_${imag}.txt"
		./vfem_pml | tee "full_${real}_${imag}.txt"
		./vfem_pml_small | tee "small_${real}_${imag}.txt"
		pltname=`ls loop_pml*.plt | sort | tail -1`
		pltname_new=`echo ${pltname} | sed "s/\./\.${real}_${imag}\./g"`
		mv ${pltname} ${pltname_new}
		datname=`ls loop_pml*.dat | sort | tail -1`
		datname_new=`echo ${datname} | sed "s/\./\.${real}_${imag}\./g"`
		mv ${datname} ${datname_new}
	done
done

