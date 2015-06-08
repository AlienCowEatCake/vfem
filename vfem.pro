TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += \
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

HEADERS += \
    src/common/matrix.h \
    src/common/common.h \
    src/geometry/point.h \
    src/geometry/vector3.h \
    src/elements/edge.h \
    src/elements/face.h \
    src/elements/triangle.h \
    src/elements/tetrahedron.h \
    src/common/basis_config.h \
    src/common/integration.h \
    src/elements/octal_tree.h \
    src/vfem/phys.h \
    src/vfem/slae.h \
    src/vfem/vfem.h \
    src/problems/problems.h \
    src/solvers/COCG_LLT_Smooth.h

CONFIG += warn_on

*g++*|*clang* {
    QMAKE_CXXFLAGS *= -ansi
    QMAKE_CXXFLAGS *= -pedantic
    QMAKE_CXXFLAGS_WARN_ON *= -Wextra
    QMAKE_CXXFLAGS_RELEASE -= -O2
    QMAKE_CXXFLAGS_RELEASE *= -O3
    QMAKE_CXXFLAGS_RELEASE *= -march=native
    QMAKE_CXXFLAGS_RELEASE *= -mtune=native
    QMAKE_CXXFLAGS_RELEASE *= -DNDEBUG
}

*msvc* {
    QMAKE_CXXFLAGS_RELEASE -= -O2
    QMAKE_CXXFLAGS_RELEASE *= -Ox
    DEFINES += _CRT_SECURE_NO_WARNINGS
    DEFINES += _USE_MATH_DEFINES
}
