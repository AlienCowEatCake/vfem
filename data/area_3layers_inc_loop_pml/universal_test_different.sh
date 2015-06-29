#!/bin/bash

INP_FILE="template_mesh3_inc_z=-5.geo"
OUT_FILE="autogen_mesh3_inc_z=-5.geo"
OUT_FILE_1="autogen_mesh3_inc_z=-5_small.msh"
OUT_FILE_2="autogen_mesh3_inc_z=-5_full.msh"
MESH_DIR="data/area_3layers_inc_loop_pml"
CURR_DIR=`pwd`
GMSH="gmsh"

function calc {

cd "${CURR_DIR}"
mkdir -p "${RESULT_DIR}/${MESH_DIR}"
cp -a "${CURR_DIR}/${MESH_DIR}/*.geo" "${CURR_DIR}/${MESH_DIR}/*.txt" "${CURR_DIR}/${RESULT_DIR}/${MESH_DIR}/"
cd "${CURR_DIR}/${RESULT_DIR}/${MESH_DIR}"

cat "${INP_FILE}" | sed "s/PML_BEGIN/${PML_BEGIN}/g ; s/PML_WIDTH/${PML_WIDTH}/g" |
sed "s/FIRST_PHASE_HERE/Mesh.Optimize = 1;\nMesh.OptimizeNetgen = 1;\nMesh 3;\nSave \"${OUT_FILE_1}\";\n/g" |
sed "s/SECOND_PHASE_HERE/Mesh 1;\nMesh.Optimize = 1;\nMesh.OptimizeNetgen = 1;\nMesh 3;\nSave \"${OUT_FILE_2}\";\n/g" > "${OUT_FILE}"
"${GMSH}" - "${OUT_FILE}"
if [[ $? -ne 0 ]] ; then
	cat "${INP_FILE}" | sed "s/PML_BEGIN/${PML_BEGIN}/g ; s/PML_WIDTH/${PML_WIDTH}/g" |
	sed "s/FIRST_PHASE_HERE/Mesh.Optimize = 1;\nMesh.OptimizeNetgen = 0;\nMesh 3;\nSave \"${OUT_FILE_1}\";\n/g" |
	sed "s/SECOND_PHASE_HERE/Mesh 1;\nMesh.Optimize = 1;\nMesh.OptimizeNetgen = 0;\nMesh 3;\nSave \"${OUT_FILE_2}\";\n/g" > "${OUT_FILE}"
	"${GMSH}" - "${OUT_FILE}"
	if [[ $? -ne 0 ]] ; then
		cat "${INP_FILE}" | sed "s/PML_BEGIN/${PML_BEGIN}/g ; s/PML_WIDTH/${PML_WIDTH}/g" |
		sed "s/FIRST_PHASE_HERE/Mesh.Optimize = 0;\nMesh.OptimizeNetgen = 0;\nMesh 3;\nSave \"${OUT_FILE_1}\";\n/g" |
		sed "s/SECOND_PHASE_HERE/Mesh 1;\nMesh.Optimize = 0;\nMesh.OptimizeNetgen = 0;\nMesh 3;\nSave \"${OUT_FILE_2}\";\n/g" > "${OUT_FILE}"
		"${GMSH}" - "${OUT_FILE}"
	fi
fi

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

exit 0
