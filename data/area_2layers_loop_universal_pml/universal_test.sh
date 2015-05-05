#!/bin/bash

PML_BEGIN=600
PML_WIDTH=100
PML_CHI_RE=4
PML_CHI_IM=1
PML_M=3
INP_FILE="autogen2_template.geo"
OUT_FILE="autogen2.geo"
OUT_FILE_1="universal_small.msh"
OUT_FILE_2="universal_full.msh"
MESH_DIR="data/area_2layers_loop_universal_pml"
CURR_DIR=`pwd`
GMSH="gmsh"
RESULT_DIR="test_begin=${PML_BEGIN}_width=${PML_WIDTH}_chi=(${PML_CHI_RE},${PML_CHI_IM})_m=${PML_M}"

#set -x

cd "${MESH_DIR}"
rm -f "${OUT_FILE_1}" "${OUT_FILE_2}" 2>/dev/null

cat "${INP_FILE}" | sed "s/PML_BEGIN/${PML_BEGIN}/g ; s/PML_WIDTH/${PML_WIDTH}/g" |
sed "s/FIRST_PHASE_HERE/Mesh.Optimize = 1;\nMesh.OptimizeNetgen = 1;\nMesh 3;\nSave Sprintf(\"${OUT_FILE_1}\", 1);\n/g" |
sed "s/SECOND_PHASE_HERE/Mesh.Optimize = 1;\nMesh.OptimizeNetgen = 1;\nMesh 3;\nSave Sprintf(\"${OUT_FILE_2}\", 1);\n/g" > "${OUT_FILE}"
"${GMSH}" - "${OUT_FILE}"
if [[ $? -ne 0 ]] ; then
	cat "${INP_FILE}" | sed "s/PML_BEGIN/${PML_BEGIN}/g ; s/PML_WIDTH/${PML_WIDTH}/g" |
	sed "s/FIRST_PHASE_HERE/Mesh.Optimize = 1;\nMesh.OptimizeNetgen = 0;\nMesh 3;\nSave Sprintf(\"${OUT_FILE_1}\", 1);\n/g" |
	sed "s/SECOND_PHASE_HERE/Mesh.Optimize = 1;\nMesh.OptimizeNetgen = 0;\nMesh 3;\nSave Sprintf(\"${OUT_FILE_2}\", 1);\n/g" > "${OUT_FILE}"
	"${GMSH}" - "${OUT_FILE}"
	if [[ $? -ne 0 ]] ; then
		cat "${INP_FILE}" | sed "s/PML_BEGIN/${PML_BEGIN}/g ; s/PML_WIDTH/${PML_WIDTH}/g" |
		sed "s/FIRST_PHASE_HERE/Mesh.Optimize = 0;\nMesh.OptimizeNetgen = 0;\nMesh 3;\nSave Sprintf(\"${OUT_FILE_1}\", 1);\n/g" |
		sed "s/SECOND_PHASE_HERE/Mesh.Optimize = 0;\nMesh.OptimizeNetgen = 0;\nMesh 3;\nSave Sprintf(\"${OUT_FILE_2}\", 1);\n/g" > "${OUT_FILE}"
		"${GMSH}" - "${OUT_FILE}"
	fi
fi

cd "${CURR_DIR}"
mkdir -p "${RESULT_DIR}"
cd "${RESULT_DIR}"
ln -s "${CURR_DIR}/data" "data"
echo -e "${PML_BEGIN}\n\n" > "pml_begin.txt"
echo -e "${PML_WIDTH}\n\n" > "pml_width.txt"
echo -e "${PML_CHI_RE} ${PML_CHI_IM}\n\n" > "pml_chi.txt"
echo -e "${PML_M}\n\n" > "pml_m.txt"
"${CURR_DIR}/vfem_classic" | tee "out_classic.txt"
"${CURR_DIR}/vfem_pml" | tee "out_pml.txt"
"${CURR_DIR}/vfem_pml_small" | tee "out_pml_small.txt"

exit 0
