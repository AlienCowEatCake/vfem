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
	src/vfem/slae.cpp \
	src/vfem/vfem.cpp \
	src/vfem/vfem_input.cpp \
	src/vfem/vfem_output.cpp \
	src/problems/analytical_cube.cpp \
	src/problems/field_cube.cpp \
	src/solvers/COCG_LLT_Smooth.cpp
OBJECTS = $(SOURCES:.cpp=.o)

all: $(SOURCES) $(EXECUTABLE)

.PHONY: clean install

$(EXECUTABLE): $(OBJECTS) 
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJECTS) $(EXECUTABLE)

