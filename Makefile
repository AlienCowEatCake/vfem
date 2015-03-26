CXX ?= g++
CXXFLAGS_EXTRA ?= -Wall -Wextra -ansi -pedantic -pipe
CXXFLAGS_OPTIMIZE ?= -O3 -march=native -mtune=native -DNDEBUG
LDFLAGS_EXTRA ?= -s
EXECUTABLE = vfem

LINK.o = $(LINK.cc)
CXXFLAGS += $(CXXFLAGS_EXTRA) $(CXXFLAGS_OPTIMIZE)
LDFLAGS += $(LDFLAGS_EXTRA)

SOURCES = \
	src/main.cpp \
	src/geometry/point.cpp \
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
	src/problems/analytical_cube.cpp \
	src/problems/area_pml_source.cpp \
	src/problems/source_pml.cpp \
	src/solvers/BiCGComplex_VC.cpp \
	src/solvers/BiCGStabComplex_VC.cpp \
	src/solvers/CGMComplex_LLT.cpp \
	src/solvers/CGMComplex_VC.cpp
OBJECTS = $(SOURCES:.cpp=.o)

all: $(SOURCES) $(EXECUTABLE)

.PHONY: clean install

$(EXECUTABLE): $(OBJECTS) 
	$(CXX) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJECTS) $(EXECUTABLE)

