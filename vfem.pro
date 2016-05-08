# Пример генерации проекта для MSVC:
# set QMAKESPEC=win32-msvc2008
# qmake -tp vc
#
# Пример получения списка cpp-шников для заполнения Makefile:
# cat src/core/core.pri vfem.pro | grep '.cpp' | sed 's/$$PWD/src\/core/g ; s/^    /\t/g ; s/cpp$/cpp \\/'
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
    src/vfem/config/config.cpp \
    src/vfem/common/cubatures.cpp \
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
    src/vfem/config/config.h \
    src/vfem/config/evaluator_helmholtz.h \
    src/vfem/common/common.h \
    src/vfem/common/cubatures.h \
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
