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
    src/elements/tetrahedron_pml.cpp \
    src/vfem/slae.cpp \
    src/vfem/vfem.cpp \
    src/vfem/vfem_input.cpp \
    src/vfem/vfem_output.cpp \
    src/problems/analytical_cube.cpp \
    src/problems/loop_pml.cpp \
    src/problems/area_2layers_loop_pml.cpp \
#    src/solvers/BiCGComplex_VC.cpp \
#    src/solvers/BiCGStabComplex_VC.cpp \
#    src/solvers/CGMComplex_LLT.cpp \
#    src/solvers/CGMComplex_VC.cpp \
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
    src/elements/octal_tree.h \
    src/vfem/phys.h \
    src/vfem/slae.h \
    src/vfem/vfem.h \
    src/problems/problems.h \
#    src/solvers/BiCGComplex_VC.h \
#    src/solvers/BiCGStabComplex_VC.h \
#    src/solvers/CGMComplex_LLT.h \
#    src/solvers/CGMComplex_VC.h \
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
}

*msvc* {
    QMAKE_CXXFLAGS_RELEASE -= -O2
    QMAKE_CXXFLAGS_RELEASE *= -Ox
    DEFINES += _CRT_SECURE_NO_WARNINGS
    DEFINES += _USE_MATH_DEFINES
    # Добавлено: 24 Feb 2015
    # MSVC 2013 использует странные оптимизации
    *msvc2013* {
        QMAKE_CXXFLAGS_RELEASE -= -Ox
        QMAKE_CXXFLAGS_RELEASE *= -O1
    }
}
