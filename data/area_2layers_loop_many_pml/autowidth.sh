#!/bin/bash

for i in `LANG=C seq 50 10 200`
do
	echo "==============================================================================="
	cd "data/area_2layers_loop_many_pml"
	position=$((${i}+600))
	cat "autogen_template.geo" | sed "s/REPLACEME/${position}/g" > "autogen.geo"
	rm -f "autogen.msh" "area_2layers_loop_many_pml_slae.txt" 2> /dev/null
	gmsh -3 -optimize -optimize_netgen "autogen.geo"
	cd ../..
	echo -e "${i}\n" > "pml_width.txt"
	echo "==============================================================================="
	echo "Testing with width = ${i}"
	./vfem_classic | tee "classic_${i}.txt"
	./vfem_pml | tee "pml_${i}.txt"
done

for i in `ls pml_*.txt | grep -v width`
do
	echo -n ${i} | sed 's/.*_//g ; s/\..*//g'
	echo -n " "
	cat ${i} | grep Diff | sed 's/.* //g'
done | sort --human-numeric-sort | tee width.txt
