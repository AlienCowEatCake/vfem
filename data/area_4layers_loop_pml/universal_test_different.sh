#!/bin/bash

INP_FILE="mesh4_z=-5.geo"
OUT_FILE_1="mesh4_z=-5_small.msh"
OUT_FILE_2="mesh4_z=-5_full.msh"
MESH_DIR="data/area_4layers_loop_pml"
CURR_DIR=`pwd`
GMSH="nosrand gmsh"
SLAE_FILE="area_4layers_loop_pml_slae.txt"
GOOD_DIFF="0.1"

# ${GMSH} -3 -optimize -o "${OUT_FILE_2}" "${INP_FILE}"
# "${CURR_DIR}/"crop_mesh "${OUT_FILE_2}" "${OUT_FILE_1}" 31 32 33 34 52
# "${CURR_DIR}/"vfem_classic | tee "out_classic.txt"

##less.c
##include <stdio.h>
#int main(int argc, char * argv[]) {
#double a, b;
#if(argc != 3) return 1;
#sscanf(argv[1], "%lf", &a);
#sscanf(argv[2], "%lf", &b);
#printf("%d\n", (a < b ? 1 : 0));
#return 0; }


function calc_nomesh_notime {

cd "${CURR_DIR}"
#if [ -d "${RESULT_DIR}/${MESH_DIR}" ]; then return; fi
mkdir -p "${RESULT_DIR}/${MESH_DIR}"
cp -a "${CURR_DIR}/${MESH_DIR}/"*.msh "${CURR_DIR}/${MESH_DIR}/"*.txt "${CURR_DIR}/${RESULT_DIR}/${MESH_DIR}/"
cp -a "${CURR_DIR}/${SLAE_FILE}" "${CURR_DIR}/${RESULT_DIR}/"
cd "${CURR_DIR}/${RESULT_DIR}"
echo -e "${PML_BEGIN}\n\n" > "pml_begin.txt"
echo -e "${PML_WIDTH}\n\n" > "pml_width.txt"
echo -e "${PML_CHI_AIR_RE} ${PML_CHI_AIR_IM}\n\n" > "pml_chi_air.txt"
echo -e "${PML_CHI_LAYER1_RE} ${PML_CHI_LAYER1_IM}\n\n" > "pml_chi_layer1.txt"
echo -e "${PML_CHI_LAYER2_RE} ${PML_CHI_LAYER2_IM}\n\n" > "pml_chi_layer2.txt"
echo -e "${PML_CHI_LAYER3_RE} ${PML_CHI_LAYER3_IM}\n\n" > "pml_chi_layer3.txt"
echo -e "${PML_M}\n\n" > "pml_m.txt"
"${CURR_DIR}/vfem_pml" | tee "out_pml.txt"
diff=`cat "out_pml.txt" | grep 'Diff (L2)' | sed 's/.*\t//g'`
if [ `"${CURR_DIR}/less" "${diff}" "${GOOD_DIFF}"` -eq 1 ] ; then
"${CURR_DIR}/vfem_pml_small" | tee "out_pml_small.txt"
fi
cd "${CURR_DIR}"

}

function calc_nomesh {

cd "${CURR_DIR}"
#if [ -d "${RESULT_DIR}/${MESH_DIR}" ]; then return; fi
mkdir -p "${RESULT_DIR}/${MESH_DIR}"
cp -a "${CURR_DIR}/${MESH_DIR}/"*.msh "${CURR_DIR}/${MESH_DIR}/"*.txt "${CURR_DIR}/${RESULT_DIR}/${MESH_DIR}/"
cp -a "${CURR_DIR}/${SLAE_FILE}" "${CURR_DIR}/${RESULT_DIR}/"
cd "${CURR_DIR}/${RESULT_DIR}"
echo -e "${PML_BEGIN}\n\n" > "pml_begin.txt"
echo -e "${PML_WIDTH}\n\n" > "pml_width.txt"
echo -e "${PML_CHI_AIR_RE} ${PML_CHI_AIR_IM}\n\n" > "pml_chi_air.txt"
echo -e "${PML_CHI_LAYER1_RE} ${PML_CHI_LAYER1_IM}\n\n" > "pml_chi_layer1.txt"
echo -e "${PML_CHI_LAYER2_RE} ${PML_CHI_LAYER2_IM}\n\n" > "pml_chi_layer2.txt"
echo -e "${PML_CHI_LAYER3_RE} ${PML_CHI_LAYER3_IM}\n\n" > "pml_chi_layer3.txt"
echo -e "${PML_M}\n\n" > "pml_m.txt"
"${CURR_DIR}/vfem_pml" | tee "out_pml.txt"
"${CURR_DIR}/vfem_pml_small" | tee "out_pml_small.txt"
cd "${CURR_DIR}"

}

PML_BEGIN=600
PML_WIDTH=100
PML_M=3
for j in `LANG=C seq 0 1`
do
for i in `LANG=C seq 2 6`
do
PML_CHI_AIR_RE=$i
PML_CHI_AIR_IM=$j
PML_CHI_LAYER1_RE=$i
PML_CHI_LAYER1_IM=$j
PML_CHI_LAYER2_RE=$i
PML_CHI_LAYER2_IM=$j
PML_CHI_LAYER3_RE=$i
PML_CHI_LAYER3_IM=$j
RESULT_DIR="test_begin=${PML_BEGIN}_width=${PML_WIDTH}_chi-air=(${PML_CHI_AIR_RE},${PML_CHI_AIR_IM})_chi-l1=(${PML_CHI_LAYER1_RE},${PML_CHI_LAYER1_IM})_chi-l2=(${PML_CHI_LAYER2_RE},${PML_CHI_LAYER2_IM})_chi-l3=(${PML_CHI_LAYER3_RE},${PML_CHI_LAYER3_IM})_m=${PML_M}"
calc_nomesh_notime
done
done
exit 0

PML_BEGIN=600
PML_WIDTH=100
PML_CHI_AIR_RE=5
PML_CHI_AIR_IM=1
PML_CHI_LAYER1_RE=2
PML_CHI_LAYER1_IM=2
PML_CHI_LAYER2_RE=2
PML_CHI_LAYER2_IM=2
PML_CHI_LAYER3_RE=2
PML_CHI_LAYER3_IM=2
PML_M=3
RESULT_DIR="test_begin=${PML_BEGIN}_width=${PML_WIDTH}_chi-air=(${PML_CHI_AIR_RE},${PML_CHI_AIR_IM})_chi-l1=(${PML_CHI_LAYER1_RE},${PML_CHI_LAYER1_IM})_chi-l2=(${PML_CHI_LAYER2_RE},${PML_CHI_LAYER2_IM})_chi-l3=(${PML_CHI_LAYER3_RE},${PML_CHI_LAYER3_IM})_m=${PML_M}"
#calc_nomesh_notime

for PML_CHI_AIR_RE in `LANG=C seq 3 6` ; do
for PML_CHI_AIR_IM in `LANG=C seq 0 1` ; do
for PML_CHI_LAYER1_RE in `LANG=C seq 0 1` ; do
for PML_CHI_LAYER1_IM in `LANG=C seq 5 7` ; do
for PML_CHI_LAYER2_RE in `LANG=C seq 1 3` ; do
for PML_CHI_LAYER2_IM in `LANG=C seq 1 3` ; do
for PML_CHI_LAYER3_RE in `LANG=C seq 1 3` ; do
for PML_CHI_LAYER3_IM in `LANG=C seq 1 3` ; do
RESULT_DIR="test_begin=${PML_BEGIN}_width=${PML_WIDTH}_chi-air=(${PML_CHI_AIR_RE},${PML_CHI_AIR_IM})_chi-l1=(${PML_CHI_LAYER1_RE},${PML_CHI_LAYER1_IM})_chi-l2=(${PML_CHI_LAYER2_RE},${PML_CHI_LAYER2_IM})_chi-l3=(${PML_CHI_LAYER3_RE},${PML_CHI_LAYER3_IM})_m=${PML_M}"
#calc_nomesh_notime
done
done
done
done
done
done
done
done

exit 0

