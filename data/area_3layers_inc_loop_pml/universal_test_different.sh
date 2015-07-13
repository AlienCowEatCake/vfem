#!/bin/bash

INP_FILE="template_mesh3_inc_z=-5.geo"
OUT_FILE="autogen_mesh3_inc_z=-5.geo"
OUT_FILE_1="autogen_mesh3_inc_z=-5_small.msh"
OUT_FILE_2="autogen_mesh3_inc_z=-5_full.msh"
MESH_DIR="data/area_3layers_inc_loop_pml"
CURR_DIR=`pwd`
GMSH="nosrand gmsh"
SLAE_FILE="area_3layers_inc_loop_pml_slae.txt"

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
echo -e "${PML_CHI_WATER_RE} ${PML_CHI_WATER_IM}\n\n" > "pml_chi_water.txt"
echo -e "${PML_CHI_GROUND_RE} ${PML_CHI_GROUND_IM}\n\n" > "pml_chi_ground.txt"
echo -e "${PML_M}\n\n" > "pml_m.txt"
"${CURR_DIR}/vfem_pml" | tee "out_pml.txt"
"${CURR_DIR}/vfem_pml_small" | tee "out_pml_small.txt"
cd "${CURR_DIR}"

}

function calc {

cd "${CURR_DIR}"
#if [ -d "${RESULT_DIR}/${MESH_DIR}" ]; then return; fi
mkdir -p "${RESULT_DIR}/${MESH_DIR}"
cp -a "${CURR_DIR}/${MESH_DIR}/"*.geo "${CURR_DIR}/${MESH_DIR}/"*.txt "${CURR_DIR}/${RESULT_DIR}/${MESH_DIR}/"
cd "${CURR_DIR}/${RESULT_DIR}/${MESH_DIR}"

#cat "${INP_FILE}" | sed "s/PML_BEGIN/${PML_BEGIN}/g ; s/PML_WIDTH/${PML_WIDTH}/g" |
#sed "s/FIRST_PHASE_HERE/Mesh.Optimize = 1;\nMesh.OptimizeNetgen = 1;\nMesh 3;\nSave \"${OUT_FILE_1}\";\n/g" |
#sed "s/SECOND_PHASE_HERE/Mesh 1;\nMesh.Optimize = 1;\nMesh.OptimizeNetgen = 1;\nMesh 3;\nSave \"${OUT_FILE_2}\";\n/g" > "${OUT_FILE}"
#${GMSH} - "${OUT_FILE}"
cat "${INP_FILE}" | sed "s/PML_BEGIN/${PML_BEGIN}/g ; s/PML_WIDTH/${PML_WIDTH}/g" |
sed "s/FIRST_PHASE_HERE//g ; s/SECOND_PHASE_HERE//g" > "${OUT_FILE}"
${GMSH} -3 -optimize -optimize_netgen -o "${OUT_FILE_2}" "${OUT_FILE}"
if [[ $? -ne 0 ]] ; then
	#cat "${INP_FILE}" | sed "s/PML_BEGIN/${PML_BEGIN}/g ; s/PML_WIDTH/${PML_WIDTH}/g" |
	#sed "s/FIRST_PHASE_HERE/Mesh.Optimize = 1;\nMesh.OptimizeNetgen = 0;\nMesh 3;\nSave \"${OUT_FILE_1}\";\n/g" |
	#sed "s/SECOND_PHASE_HERE/Mesh 1;\nMesh.Optimize = 1;\nMesh.OptimizeNetgen = 0;\nMesh 3;\nSave \"${OUT_FILE_2}\";\n/g" > "${OUT_FILE}"
	#${GMSH} - "${OUT_FILE}"
	${GMSH} -3 -optimize -o "${OUT_FILE_2}" "${OUT_FILE}"
	if [[ $? -ne 0 ]] ; then
		#cat "${INP_FILE}" | sed "s/PML_BEGIN/${PML_BEGIN}/g ; s/PML_WIDTH/${PML_WIDTH}/g" |
		#sed "s/FIRST_PHASE_HERE/Mesh.Optimize = 0;\nMesh.OptimizeNetgen = 0;\nMesh 3;\nSave \"${OUT_FILE_1}\";\n/g" |
		#sed "s/SECOND_PHASE_HERE/Mesh 1;\nMesh.Optimize = 0;\nMesh.OptimizeNetgen = 0;\nMesh 3;\nSave \"${OUT_FILE_2}\";\n/g" > "${OUT_FILE}"
		#${GMSH} - "${OUT_FILE}"
		${GMSH} -3 -o "${OUT_FILE_2}" "${OUT_FILE}"
	fi
fi
"${CURR_DIR}/crop_mesh" "${OUT_FILE_2}" "${OUT_FILE_1}" 31 32 33 52

cd "${CURR_DIR}/${RESULT_DIR}"
echo -e "${PML_BEGIN}\n\n" > "pml_begin.txt"
echo -e "${PML_WIDTH}\n\n" > "pml_width.txt"
echo -e "${PML_CHI_AIR_RE} ${PML_CHI_AIR_IM}\n\n" > "pml_chi_air.txt"
echo -e "${PML_CHI_WATER_RE} ${PML_CHI_WATER_IM}\n\n" > "pml_chi_water.txt"
echo -e "${PML_CHI_GROUND_RE} ${PML_CHI_GROUND_IM}\n\n" > "pml_chi_ground.txt"
echo -e "${PML_M}\n\n" > "pml_m.txt"
"${CURR_DIR}/vfem_classic" | tee "out_classic.txt"
"${CURR_DIR}/vfem_pml" | tee "out_pml.txt"
"${CURR_DIR}/vfem_pml_small" | tee "out_pml_small.txt"

cd "${CURR_DIR}"

}

PML_BEGIN=600
PML_WIDTH=100
PML_CHI_AIR_RE=5
PML_CHI_AIR_IM=1
PML_CHI_WATER_RE=0
PML_CHI_WATER_IM=7
PML_CHI_GROUND_RE=2
PML_CHI_GROUND_IM=2
PML_M=3
RESULT_DIR="test_begin=${PML_BEGIN}_width=${PML_WIDTH}_chi-air=(${PML_CHI_AIR_RE},${PML_CHI_AIR_IM})_chi-water=(${PML_CHI_WATER_RE},${PML_CHI_WATER_IM})_chi-ground=(${PML_CHI_GROUND_RE},${PML_CHI_GROUND_IM})_m=${PML_M}"
calc

for PML_CHI_AIR_RE in `LANG=C seq 3 6` ; do
for PML_CHI_AIR_IM in `LANG=C seq 0 1` ; do
for PML_CHI_WATER_RE in `LANG=C seq 0 1` ; do
for PML_CHI_WATER_IM in `LANG=C seq 5 7` ; do
for PML_CHI_GROUND_RE in `LANG=C seq 1 3` ; do
for PML_CHI_GROUND_IM in `LANG=C seq 1 3` ; do
RESULT_DIR="test_begin=${PML_BEGIN}_width=${PML_WIDTH}_chi-air=(${PML_CHI_AIR_RE},${PML_CHI_AIR_IM})_chi-water=(${PML_CHI_WATER_RE},${PML_CHI_WATER_IM})_chi-ground=(${PML_CHI_GROUND_RE},${PML_CHI_GROUND_IM})_m=${PML_M}"
calc
done
done
done
done
done
done

exit 0
