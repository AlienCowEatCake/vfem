CXX ?= g++
CXXFLAGS_EXTRA ?= -Wall -Wextra -std=c++0x -pedantic -pipe -DUSE_NOSIGHUP
CXXFLAGS_OPTIMIZE ?= -O3 -march=native -mtune=native -DNDEBUG
LDFLAGS_EXTRA ?= -Wl,-O1 -s
EXECUTABLE = vfem

LINK.o = $(LINK.cc)
CXXFLAGS += $(CXXFLAGS_EXTRA) $(CXXFLAGS_OPTIMIZE) -fopenmp -DUSE_OMP
LDFLAGS += $(LDFLAGS_EXTRA) -lrt -fopenmp

SOURCES = \
	src/core/cubatures/line_integration.cpp \
	src/core/cubatures/tetrahedron_integration.cpp \
	src/core/cubatures/triangle_integration.cpp \
	src/core/evaluator/evaluator_internal/jit/common.cpp \
	src/core/evaluator/evaluator_internal/jit/func_templates.cpp \
	src/core/evaluator/evaluator_internal/jit/oper_templates.cpp \
	src/core/evaluator/evaluator_internal/jit/real_templates.cpp \
	src/core/evaluator/evaluator_internal/transition_table.cpp \
	src/core/solvers/CSLR/solvers_factory.cpp \
	src/core/solvers/CSLR/symmetric/complex/BiCG_Complex/BiCG_Complex.cpp \
	src/core/solvers/CSLR/symmetric/complex/BiCG_Complex/BiCG_Complex_Smooth.cpp \
	src/core/solvers/CSLR/symmetric/complex/BiCGStab_Complex/BiCGStab_Complex.cpp \
	src/core/solvers/CSLR/symmetric/complex/BiCGStab_Complex/BiCGStab_Complex_Smooth.cpp \
	src/core/solvers/CSLR/symmetric/complex/COCG/COCG.cpp \
	src/core/solvers/CSLR/symmetric/complex/COCG/COCG_Smooth.cpp \
	src/core/solvers/CSLR/symmetric/complex/COCG/COCG_Smooth_MKL.cpp \
	src/core/solvers/CSLR/symmetric/complex/COCG/COCG_Smooth_OpenMP.cpp \
	src/core/solvers/CSLR/symmetric/complex/COCR/COCR.cpp \
	src/core/solvers/CSLR/symmetric/complex/COCR/COCR_Smooth.cpp \
	src/core/solvers/CSLR/symmetric/complex/GMRES_Complex/GMRES_Complex.cpp \
	src/core/solvers/CSLR/symmetric/complex/GMRES_Complex/GMRES_Complex_MKL.cpp \
	src/core/solvers/CSLR/symmetric/complex/GMRES_Complex/GMRES_Complex_OpenMP.cpp \
	src/core/utils/inifile.cpp \
	src/core/utils/nosighup.cpp \
	src/core/utils/progress.cpp \
	src/core/utils/strings.cpp \
	src/core/utils/timers.cpp \
	src/core/wrappers/mkl_wrapper.cpp \
	src/core/wrappers/omp_wrapper.cpp \
	src/vfem/common/config.cpp \
	src/vfem/elements/edge.cpp \
	src/vfem/elements/tetrahedron.cpp \
	src/vfem/elements/tetrahedron_pml.cpp \
	src/vfem/elements/triangle.cpp \
	src/vfem/main.cpp \
	src/vfem/problems/standard.cpp \
	src/vfem/problems/standard_diff.cpp \
	src/vfem/problems/standard_pml.cpp \
	src/vfem/vfem/phys.cpp \
	src/vfem/vfem/slae.cpp \
	src/vfem/vfem/vfem.cpp \
	src/vfem/vfem/vfem_input.cpp \
	src/vfem/vfem/vfem_output.cpp \
	src/vfem/vfem/vfem_v_cycle.cpp
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

