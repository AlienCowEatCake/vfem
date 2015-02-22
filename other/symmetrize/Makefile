CXX ?= g++
CXXFLAGS_EXTRA ?= -Wall -Wextra -ansi -pedantic
CXXFLAGS_OPTIMIZE ?= -O3 -march=native -mtune=native
LDFLAGS_EXTRA ?= -s
EXECUTABLE = symmetrize

LINK.o = $(LINK.cc)
CXXFLAGS += $(CXXFLAGS_EXTRA) $(CXXFLAGS_OPTIMIZE)
LDFLAGS += $(LDFLAGS_EXTRA)

SOURCES = \
	symmetrize.cpp
OBJECTS = $(SOURCES:.cpp=.o)

all: $(SOURCES) $(EXECUTABLE)

.PHONY: clean install

$(EXECUTABLE): $(OBJECTS) 
	$(CXX) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJECTS) $(EXECUTABLE)
