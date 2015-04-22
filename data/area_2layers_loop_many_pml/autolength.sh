#!/bin/bash

for i in `LANG=C seq 500 10 1400`
do
	echo "==============================================================================="
	cd "data/area_2layers_loop_many_pml"
	cat "autogen_template.geo" | sed "s/REPLACEME/${i}/g" > "autogen2.geo"
	rm -f "autogen2.msh" "area_2layers_loop_many_pml_slae.txt" 2> /dev/null
	gmsh -3 -optimize -optimize_netgen "autogen2.geo"
	cd ../..
	echo "==============================================================================="
	echo "Testing with length = ${i}"
	./vfem_classic | tee "classic_${i}.txt"
	./vfem_pml | tee "pml_${i}.txt"
done

for i in `ls pml_*.txt`
do
	echo -n ${i} | sed 's/.*_//g ; s/\..*//g'
	echo -n " "
	cat ${i} | grep Diff | sed 's/.* //g'
done | sort --human-numeric-sort | tee length.txt
