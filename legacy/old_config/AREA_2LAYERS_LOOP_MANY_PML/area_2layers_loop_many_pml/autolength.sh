#!/bin/bash

out_geofile="autogen2.geo"
for i in `LANG=C seq 500 100 1400`
do
	echo "==============================================================================="
	cd "data/area_2layers_loop_many_pml"
	cat "autogen2_template.geo" | sed "s/REPLACEME/${i}/g" > "${out_geofile}"
	rm -f "autogen2.msh" "area_2layers_loop_many_pml_slae.txt" 2> /dev/null
	gmsh -3 -optimize -optimize_netgen "${out_geofile}"
	if [[ $? -ne 0 ]] ; then
		gmsh -3 -optimize "${out_geofile}"
		if [[ $? -ne 0 ]] ; then
			gmsh -3 "${out_geofile}"
		fi
	fi
	cd ../..
	echo "==============================================================================="
	echo "Testing with length = ${i}"
	./vfem_classic | tee "classic_${i}.txt"
	./vfem_pml | tee "pml_${i}.txt"
	mv "line_std.txt" "line_std_${i}.txt"
	mv "line_pml.txt" "line_pml_${i}.txt"
done

for i in `ls pml_*.txt`
do
	echo -n ${i} | sed 's/.*_//g ; s/\..*//g'
	echo -n " "
	cat ${i} | grep Diff | sed 's/.* //g'
done | sort --human-numeric-sort | tee length.txt
