CXX ?= g++
CXXFLAGS_EXTRA ?= -Wall -Wextra -std=c++0x -pedantic -pipe -DUSE_NOSIGHUP
CXXFLAGS_OPTIMIZE ?= -O3 -march=native -mtune=native -DNDEBUG
LDFLAGS_EXTRA ?= -s
EXECUTABLE = vfem

LINK.o = $(LINK.cc)
CXXFLAGS += $(CXXFLAGS_EXTRA) $(CXXFLAGS_OPTIMIZE) -fopenmp
LDFLAGS += $(LDFLAGS_EXTRA) -lrt

SOURCES = \
	src/main.cpp \
	src/config/config.cpp \
	src/config/evaluator/evaluator_internal/transition_table.cpp \
	src/config/evaluator/evaluator_internal/jit/common.cpp \
	src/config/evaluator/evaluator_internal/jit/func_templates.cpp \
	src/config/evaluator/evaluator_internal/jit/oper_templates.cpp \
	src/config/evaluator/evaluator_internal/jit/real_templates.cpp \
	src/common/cubatures.cpp \
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
	src/solvers/COCG_LLT_Smooth.cpp \
	src/solvers/COCG_LLT_Smooth_MKL.cpp
OBJECTS = $(SOURCES:.cpp=.o)

all: $(SOURCES) $(EXECUTABLE)

.PHONY: clean distclean

$(EXECUTABLE): $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(OBJECTS) $(LDFLAGS) -o $@

.cpp.o:
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJECTS)

distclean: clean
	rm -f $(EXECUTABLE)

