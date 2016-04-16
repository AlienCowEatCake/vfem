CXX ?= g++
CXXFLAGS_EXTRA ?= -Wall -Wextra -std=c++0x -pedantic -pipe -DUSE_NOSIGHUP
CXXFLAGS_OPTIMIZE ?= -O3 -march=native -mtune=native -DNDEBUG
LDFLAGS_EXTRA ?= -Wl,-O1 -s
EXECUTABLE = vfem

LINK.o = $(LINK.cc)
CXXFLAGS += $(CXXFLAGS_EXTRA) $(CXXFLAGS_OPTIMIZE) -fopenmp -DUSE_OMP
LDFLAGS += $(LDFLAGS_EXTRA) -lrt -fopenmp

SOURCES = \
	src/main.cpp \
	src/config/inifile.cpp \
	src/config/config.cpp \
	src/config/evaluator/evaluator_internal/transition_table.cpp \
	src/config/evaluator/evaluator_internal/jit/common.cpp \
	src/config/evaluator/evaluator_internal/jit/func_templates.cpp \
	src/config/evaluator/evaluator_internal/jit/oper_templates.cpp \
	src/config/evaluator/evaluator_internal/jit/real_templates.cpp \
	src/common/cubatures.cpp \
	src/stubs/omp_stubs.cpp \
	src/stubs/mkl_stubs.cpp \
	src/geometry/vector3.cpp \
	src/elements/edge.cpp \
	src/elements/face.cpp \
	src/elements/triangle.cpp \
	src/elements/tetrahedron.cpp \
	src/elements/tetrahedron_pml.cpp \
	src/vfem/slae.cpp \
	src/vfem/vfem.cpp \
	src/vfem/vfem_input.cpp \
	src/vfem/vfem_output.cpp \
	src/vfem/vfem_v_cycle.cpp \
	src/problems/standard.cpp \
	src/problems/standard_diff.cpp \
	src/problems/standard_pml.cpp \
	src/solvers/BiCG_Complex/legacy/BiCGComplex_VC.cpp \
	src/solvers/BiCGStab_Complex/legacy/BiCGStabComplex_VC.cpp \
	src/solvers/COCG/legacy/CGMComplex_VC.cpp \
	src/solvers/COCG/legacy/CGMComplex_LLT.cpp \
	src/solvers/BiCG_Complex/BiCG_Complex/BiCG_Complex_Smooth.cpp \
	src/solvers/BiCGStab_Complex/BiCGStab_Complex/BiCGStab_Complex_Smooth.cpp \
	src/solvers/COCR/COCR/COCR.cpp \
	src/solvers/COCR/COCR/COCR_Di.cpp \
	src/solvers/COCR/COCR/COCR_Di_Smooth.cpp \
	src/solvers/COCR/COCR/COCR_LDLT.cpp \
	src/solvers/COCR/COCR/COCR_LDLT_Smooth.cpp \
	src/solvers/COCR/COCR/COCR_LLT.cpp \
	src/solvers/COCR/COCR/COCR_LLT_Smooth.cpp \
	src/solvers/GMRES_Complex/GMRES_Complex/GMRES_Complex.cpp \
	src/solvers/GMRES_Complex/GMRES_Complex/GMRES_Complex_Di.cpp \
	src/solvers/GMRES_Complex/GMRES_Complex/GMRES_Complex_LDLT.cpp \
	src/solvers/GMRES_Complex/GMRES_Complex/GMRES_Complex_LLT.cpp \
	src/solvers/GMRES_Complex/GMRES_Complex_OpenMP/GMRES_Complex_OpenMP.cpp \
	src/solvers/GMRES_Complex/GMRES_Complex_OpenMP/GMRES_Complex_Di_OpenMP.cpp \
	src/solvers/GMRES_Complex/GMRES_Complex_MKL/GMRES_Complex_MKL.cpp \
	src/solvers/GMRES_Complex/GMRES_Complex_MKL/GMRES_Complex_Di_MKL.cpp \
	src/solvers/GMRES_Complex/GMRES_Complex_MKL/GMRES_Complex_LLT_MKL.cpp \
	src/solvers/COCG/COCG/COCG.cpp \
	src/solvers/COCG/COCG/COCG_Di.cpp \
	src/solvers/COCG/COCG/COCG_GS.cpp \
	src/solvers/COCG/COCG/COCG_LDLT.cpp \
	src/solvers/COCG/COCG/COCG_LLT.cpp \
	src/solvers/COCG/COCG/COCG_Smooth.cpp \
	src/solvers/COCG/COCG/COCG_Di_Smooth.cpp \
	src/solvers/COCG/COCG/COCG_GS_Smooth.cpp \
	src/solvers/COCG/COCG/COCG_LLT_Smooth.cpp \
	src/solvers/COCG/COCG/COCG_LDLT_Smooth.cpp \
	src/solvers/COCG/COCG_OpenMP/COCG_Smooth_OpenMP.cpp \
	src/solvers/COCG/COCG_OpenMP/COCG_Di_Smooth_OpenMP.cpp \
	src/solvers/COCG/COCG_MKL/COCG_Smooth_MKL.cpp \
	src/solvers/COCG/COCG_MKL/COCG_Di_Smooth_MKL.cpp \
	src/solvers/COCG/COCG_MKL/COCG_LLT_Smooth_MKL.cpp
OBJECTS = $(SOURCES:.cpp=.o)

all: $(SOURCES) $(EXECUTABLE)

.PHONY: clean distclean

$(EXECUTABLE): $(OBJECTS)
	$(CXX) $(OBJECTS) $(LDFLAGS) -o $@

.cpp.o:
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJECTS)

distclean: clean
	rm -f $(EXECUTABLE)

