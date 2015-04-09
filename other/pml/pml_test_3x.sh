#!/bin/bash

for real_air in `LANG=C seq 0 1 1`
do
  for imag_air in `LANG=C seq 3 1 4`
  do
    echo -e "${real_air} ${imag_air}\n" > chi_air.txt
    for real_water in `LANG=C seq 0 1 1`
    do
      for imag_water in `LANG=C seq 6 1 7`
      do
        echo -e "${real_water} ${imag_water}\n" > chi_water.txt
        for real_ground in `LANG=C seq 0 1 7`
        do
#          for imag_ground in `LANG=C seq 0 1 7`
#          do
           imag_ground=${real_ground}
            echo -e "${real_ground} ${imag_ground}\n" > chi_ground.txt
            echo "Testing in # Air:(${real_air} , ${imag_air}), Water:(${real_water} , ${imag_water}), Ground:(${real_ground} , ${imag_ground})  ..."
            ./vfem | tee "A_${real_air}_${imag_air}_W_${real_water}_${imag_water}_G_${real_ground}_${imag_ground}.txt"
            pltname=`ls area_3layers_inc_loop_pml*.plt | sort | tail -1`
            pltname_new=`echo ${pltname} | sed "s/\./\.A_${real_air}_${imag_air}_W_${real_water}_${imag_water}_G_${real_ground}_${imag_ground}\./g"`
            mv ${pltname} ${pltname_new}
            datname=`ls area_3layers_inc_loop_pml*.dat | sort | tail -1`
            datname_new=`echo ${datname} | sed "s/\./\.A_${real_air}_${imag_air}_W_${real_water}_${imag_water}_G_${real_ground}_${imag_ground}\./g"`
            mv ${datname} ${datname_new}
#          done
        done
      done
    done
  done
done
