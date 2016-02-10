#!/bin/bash

for PML_CHI_AIR_RE in `LANG=C seq 2 4` ; do
for PML_CHI_AIR_IM in `LANG=C seq 0 1` ; do
for PML_CHI_GROUND_RE in `LANG=C seq 2 4` ; do
for PML_CHI_GROUND_IM in `LANG=C seq 1 3` ; do
RESULT_DIR="test_chi-air=(${PML_CHI_AIR_RE},${PML_CHI_AIR_IM})_chi-ground=(${PML_CHI_GROUND_RE},${PML_CHI_GROUND_IM})"
mkdir -p "${RESULT_DIR}"
cd "${RESULT_DIR}"
ln -s ../mesh2_inc_z=5_small.msh ./mesh2_inc_z=5_small.msh
ln -s ../config_pml.ini ./config_pml.ini
ln -s ../phys_pml.ini ./phys_pml.ini
ln -s ../area_2layers_inc_loop_slae.txt ./area_2layers_inc_loop_slae.txt
echo -e "[PML]\nm = 3\nwidth = 100\n" > pml.ini
echo -e "[PML.21]\nchi_real = ${PML_CHI_AIR_RE}\nchi_imag = ${PML_CHI_AIR_IM}\n" >> pml.ini
echo -e "[PML.22]\nchi_real = ${PML_CHI_GROUND_RE}\nchi_imag = ${PML_CHI_GROUND_IM}\n" >> pml.ini
../vfem -diff_simple -nosolve -nopost ../config_std.ini config_pml.ini ../diff.ini | tee out-pml.txt
cd ..
done
done
done
done

