CXX ?= g++
CXXFLAGS_EXTRA ?= -Wall -Wextra -ansi -pedantic -pipe
CXXFLAGS_OPTIMIZE ?= -O3 -DNDEBUG
LDFLAGS_EXTRA ?= -s
EXECUTABLE = crop_mesh

LINK.o = $(LINK.cc)
CXXFLAGS += $(CXXFLAGS_EXTRA) $(CXXFLAGS_OPTIMIZE)
LDFLAGS += $(LDFLAGS_EXTRA)

SOURCES = crop_mesh.cpp
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

