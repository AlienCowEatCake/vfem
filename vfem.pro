TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt
CONFIG += warn_on

SOURCES += \
    src/main.cpp \
    src/config/inifile.cpp \
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
    src/solvers/BiCGComplex_VC.cpp \
    src/solvers/BiCGStabComplex_VC.cpp \
    src/solvers/CGMComplex_LLT.cpp \
    src/solvers/CGMComplex_VC.cpp \
    src/solvers/COCG_LLT_Smooth.cpp \
    src/solvers/COCG_LLT_Smooth_MKL.cpp

HEADERS += \
    src/config/inifile.h \
    src/config/config.h \
    src/config/evaluator_helmholtz.h \
    src/config/evaluator/evaluator.h \
    src/config/evaluator/evaluator_xyz.h \
    src/config/evaluator/evaluator_operations.h \
    src/config/evaluator/evaluator_internal/type_detection.h \
    src/config/evaluator/evaluator_internal/evaluator_object.h \
    src/config/evaluator/evaluator_internal/var_container.h \
    src/config/evaluator/evaluator_internal/transition_table.h \
    src/config/evaluator/evaluator_internal/misc.h \
    src/config/evaluator/evaluator_internal/parse.h \
    src/config/evaluator/evaluator_internal/simplify.h \
    src/config/evaluator/evaluator_internal/calculate.h \
    src/config/evaluator/evaluator_internal/jit/common.h \
    src/config/evaluator/evaluator_internal/jit/opcodes.h \
    src/config/evaluator/evaluator_internal/jit/func_templates.h \
    src/config/evaluator/evaluator_internal/jit/oper_templates.h \
    src/config/evaluator/evaluator_internal/jit/real_templates.h \
    src/config/evaluator/evaluator_internal/jit/complex_templates.h \
    src/config/evaluator/evaluator_internal/jit/compile_inline.h \
    src/config/evaluator/evaluator_internal/jit/compile_extcall.h \
    src/common/trio.h \
    src/common/matrix.h \
    src/common/common.h \
    src/common/cubatures.h \
    src/geometry/point.h \
    src/geometry/vector3.h \
    src/elements/edge.h \
    src/elements/face.h \
    src/elements/triangle.h \
    src/elements/tetrahedron.h \
    src/elements/octree.h \
    src/vfem/phys.h \
    src/vfem/slae.h \
    src/vfem/vfem.h \
    src/problems/problems.h \
    src/solvers/BiCGComplex_VC.h \
    src/solvers/BiCGStabComplex_VC.h \
    src/solvers/CGMComplex_LLT.h \
    src/solvers/CGMComplex_VC.h \
    src/solvers/COCG_LLT_Smooth.h \
    src/solvers/COCG_LLT_Smooth_MKL.h \
    src/solvers/solver_interface.h

unix:!macx:QMAKE_LIBS += -lrt

*g++*|*clang* {
    QMAKE_CXXFLAGS += -fopenmp
    QMAKE_LFLAGS += -fopenmp
    QMAKE_CXXFLAGS *= -ansi
#    QMAKE_CXXFLAGS += -std=c++0x
    QMAKE_CXXFLAGS *= -pedantic
    QMAKE_CXXFLAGS_WARN_ON *= -Wextra
    QMAKE_CXXFLAGS_RELEASE -= -O2
    QMAKE_CXXFLAGS_RELEASE *= -O3
    QMAKE_CXXFLAGS_RELEASE *= -march=native
    QMAKE_CXXFLAGS_RELEASE *= -mtune=native
    QMAKE_CXXFLAGS_RELEASE *= -DNDEBUG
}

*msvc* {
    QMAKE_CXXFLAGS += -openmp
    QMAKE_CXXFLAGS_RELEASE -= -O2
    QMAKE_CXXFLAGS_RELEASE *= -Ox
    QMAKE_CXXFLAGS_RELEASE -= -GS
    QMAKE_CXXFLAGS_RELEASE *= -GS-
}

DESTDIR = .
OBJECTS_DIR = build
MOC_DIR = build
RCC_DIR = build
UI_DIR = build
