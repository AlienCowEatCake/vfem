# Пример генерации проекта для MSVC:
# set QMAKESPEC=win32-msvc2008
# qmake -tp vc
#
# Пример получения списка cpp-шников для заполнения Makefile:
# find src -name '*.cpp' | sort | sed 's/^/\t/ ; s/$/ \\/'
#
TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt
CONFIG += warn_on

CONFIG += use_mkl
CONFIG += use_omp

*g++*|*clang* {
    QMAKE_CXXFLAGS *= -ansi
    QMAKE_CXXFLAGS *= -pedantic
    QMAKE_CXXFLAGS_WARN_ON *= -Wextra
    QMAKE_CXXFLAGS_RELEASE *= -march=native
    QMAKE_CXXFLAGS_RELEASE *= -mtune=native
}

include(src/core/core.pri)

SOURCES += \
    src/vfem/main.cpp \
    src/vfem/common/config.cpp \
    src/vfem/elements/edge.cpp \
    src/vfem/elements/triangle.cpp \
    src/vfem/elements/tetrahedron.cpp \
    src/vfem/elements/tetrahedron_pml.cpp \
    src/vfem/vfem/slae.cpp \
    src/vfem/vfem/vfem.cpp \
    src/vfem/vfem/vfem_input.cpp \
    src/vfem/vfem/vfem_output.cpp \
    src/vfem/vfem/vfem_v_cycle.cpp \
    src/vfem/problems/standard.cpp \
    src/vfem/problems/standard_diff.cpp \
    src/vfem/problems/standard_pml.cpp

HEADERS += \
    src/vfem/common/config.h \
    src/vfem/common/evaluator_helmholtz.h \
    src/vfem/common/common.h \
    src/vfem/elements/edge.h \
    src/vfem/elements/triangle.h \
    src/vfem/elements/tetrahedron.h \
    src/vfem/vfem/phys.h \
    src/vfem/vfem/slae.h \
    src/vfem/vfem/vfem.h \
    src/vfem/problems/problems.h

DESTDIR = .
OBJECTS_DIR = build
MOC_DIR = build
RCC_DIR = build
UI_DIR = build
