TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt
CONFIG += warn_on

CONFIG += use_mkl
CONFIG += use_omp

unix:!macx:QMAKE_LIBS += -lrt

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
    QMAKE_CXXFLAGS_RELEASE -= -GS
    QMAKE_CXXFLAGS_RELEASE *= -GS-
}

use_omp {
    *g++* {
        QMAKE_CXXFLAGS += -fopenmp
        QMAKE_LFLAGS += -fopenmp
        DEFINES += USE_OMP
    }
    # Clang пока что не умеет OpenMP
    *msvc* {
        QMAKE_CXXFLAGS += -openmp
        DEFINES += USE_OMP
    }
}

use_mkl {
    win32-msvc201*|win32-icc* {
        contains(QMAKE_HOST.arch, x86_64) {
            PROGRAMFILES = $$system("echo %ProgramFiles(x86)%")
        } else {
            PROGRAMFILES = $$(ProgramFiles)
        }
        MKLROOT = $${PROGRAMFILES}/IntelSWTools/compilers_and_libraries_2016.2.180/windows/mkl
        exists($${MKLROOT}) {
            QMAKE_INCDIR += $$quote($${MKLROOT}/include)
            DEFINES += USE_MKL
            contains(QMAKE_TARGET.arch, x86_64) {
                QMAKE_LIBDIR += $$quote($${MKLROOT}/lib/intel64_win)
            } else {
                QMAKE_LIBDIR += $$quote($${MKLROOT}/lib/ia32_win)
            }
        }
    }
    linux-g++*-64: {
        MKLROOT = /opt/intel/compilers_and_libraries_2016.2.181/linux/mkl
        exists($${MKLROOT}) {
            MKLLIB = $${MKLROOT}/lib/intel64
            QMAKE_CXXFLAGS  -= -ansi
            QMAKE_CXXFLAGS  *= -I$${MKLROOT}/include -DUSE_MKL -std=c++0x
            QMAKE_LFLAGS    *= -Wl,-rpath -Wl,$${MKLLIB} -Wl,--no-as-needed -L$${MKLLIB}
            QMAKE_LIBS      *= -lmkl_intel_lp64 -lmkl_core -lmkl_gnu_thread
            #QMAKE_LIBS      *= -Wl,--start-group $${MKLLIB}/libmkl_intel_lp64.a $${MKLLIB}/libmkl_core.a $${MKLLIB}/libmkl_gnu_thread.a -Wl,--end-group
            QMAKE_LIBS      *= -ldl -lpthread -lm
        }
    }
}

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
    src/stubs/omp_stubs.cpp \
    src/stubs/mkl_stubs.cpp \
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
    src/solvers/BiCG_Complex/legacy/BiCGComplex_VC.cpp \
    src/solvers/BiCGStab_Complex/legacy/BiCGStabComplex_VC.cpp \
    src/solvers/COCG/legacy/CGMComplex_VC.cpp \
    src/solvers/COCG/legacy/CGMComplex_LLT.cpp \
    src/solvers/BiCG_Complex/BiCG_Complex/BiCG_Complex_Smooth.cpp \
    src/solvers/BiCGStab_Complex/BiCGStab_Complex/BiCGStab_Complex_Smooth.cpp \
    src/solvers/COCR/COCR/COCR.cpp \
    src/solvers/COCR/COCR/COCR_Di.cpp \
    src/solvers/COCR/COCR/COCR_Di_Smooth.cpp \
    src/solvers/COCR/COCR/COCR_LDLT.cpp \
    src/solvers/COCR/COCR/COCR_LDLT_Smooth.cpp \
    src/solvers/COCR/COCR/COCR_LLT.cpp \
    src/solvers/COCR/COCR/COCR_LLT_Smooth.cpp \
    src/solvers/GMRES_Complex/GMRES_Complex/GMRES_Complex.cpp \
    src/solvers/GMRES_Complex/GMRES_Complex/GMRES_Complex_Di.cpp \
    src/solvers/GMRES_Complex/GMRES_Complex/GMRES_Complex_LDLT.cpp \
    src/solvers/GMRES_Complex/GMRES_Complex/GMRES_Complex_LLT.cpp \
    src/solvers/GMRES_Complex/GMRES_Complex_OpenMP/GMRES_Complex_OpenMP.cpp \
    src/solvers/GMRES_Complex/GMRES_Complex_OpenMP/GMRES_Complex_Di_OpenMP.cpp \
    src/solvers/GMRES_Complex/GMRES_Complex_MKL/GMRES_Complex_MKL.cpp \
    src/solvers/GMRES_Complex/GMRES_Complex_MKL/GMRES_Complex_Di_MKL.cpp \
    src/solvers/GMRES_Complex/GMRES_Complex_MKL/GMRES_Complex_LLT_MKL.cpp \
    src/solvers/COCG/COCG/COCG.cpp \
    src/solvers/COCG/COCG/COCG_Di.cpp \
    src/solvers/COCG/COCG/COCG_GS.cpp \
    src/solvers/COCG/COCG/COCG_LDLT.cpp \
    src/solvers/COCG/COCG/COCG_LLT.cpp \
    src/solvers/COCG/COCG/COCG_Smooth.cpp \
    src/solvers/COCG/COCG/COCG_Di_Smooth.cpp \
    src/solvers/COCG/COCG/COCG_GS_Smooth.cpp \
    src/solvers/COCG/COCG/COCG_LLT_Smooth.cpp \
    src/solvers/COCG/COCG/COCG_LDLT_Smooth.cpp \
    src/solvers/COCG/COCG_OpenMP/COCG_Smooth_OpenMP.cpp \
    src/solvers/COCG/COCG_OpenMP/COCG_Di_Smooth_OpenMP.cpp \
    src/solvers/COCG/COCG_MKL/COCG_Smooth_MKL.cpp \
    src/solvers/COCG/COCG_MKL/COCG_Di_Smooth_MKL.cpp \
    src/solvers/COCG/COCG_MKL/COCG_LLT_Smooth_MKL.cpp

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
    src/stubs/omp_stubs.h \
    src/stubs/mkl_stubs.h \
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
    src/solvers/solver_interface.h \
    src/solvers/BiCG_Complex/legacy/BiCGComplex_VC.h \
    src/solvers/BiCGStab_Complex/legacy/BiCGStabComplex_VC.h \
    src/solvers/COCG/legacy/CGMComplex_VC.h \
    src/solvers/COCG/legacy/CGMComplex_LLT.h \
    src/solvers/BiCG_Complex/BiCG_Complex/BiCG_Complex_Smooth.h \
    src/solvers/BiCGStab_Complex/BiCGStab_Complex/BiCGStab_Complex_Smooth.h \
    src/solvers/COCR/COCR/COCR_Di.h \
    src/solvers/COCR/COCR/COCR_Di_Smooth.h \
    src/solvers/COCR/COCR/COCR.h \
    src/solvers/COCR/COCR/COCR_LDLT.h \
    src/solvers/COCR/COCR/COCR_LDLT_Smooth.h \
    src/solvers/COCR/COCR/COCR_LLT.h \
    src/solvers/COCR/COCR/COCR_LLT_Smooth.h \
    src/solvers/GMRES_Complex/GMRES_Complex/GMRES_Complex.h \
    src/solvers/GMRES_Complex/GMRES_Complex/GMRES_Complex_Di.h \
    src/solvers/GMRES_Complex/GMRES_Complex/GMRES_Complex_LDLT.h \
    src/solvers/GMRES_Complex/GMRES_Complex/GMRES_Complex_LLT.h \
    src/solvers/GMRES_Complex/GMRES_Complex_OpenMP/GMRES_Complex_OpenMP.h \
    src/solvers/GMRES_Complex/GMRES_Complex_OpenMP/GMRES_Complex_Di_OpenMP.h \
    src/solvers/GMRES_Complex/GMRES_Complex_MKL/GMRES_Complex_MKL.h \
    src/solvers/GMRES_Complex/GMRES_Complex_MKL/GMRES_Complex_Di_MKL.h \
    src/solvers/GMRES_Complex/GMRES_Complex_MKL/GMRES_Complex_LLT_MKL.h \
    src/solvers/COCG/COCG/COCG_Di.h \
    src/solvers/COCG/COCG/COCG_GS.h \
    src/solvers/COCG/COCG/COCG.h \
    src/solvers/COCG/COCG/COCG_LDLT.h \
    src/solvers/COCG/COCG/COCG_LLT.h \
    src/solvers/COCG/COCG/COCG_Smooth.h \
    src/solvers/COCG/COCG/COCG_Di_Smooth.h \
    src/solvers/COCG/COCG/COCG_GS_Smooth.h \
    src/solvers/COCG/COCG/COCG_LLT_Smooth.h \
    src/solvers/COCG/COCG/COCG_LDLT_Smooth.h \
    src/solvers/COCG/COCG_OpenMP/COCG_Smooth_OpenMP.h \
    src/solvers/COCG/COCG_OpenMP/COCG_Di_Smooth_OpenMP.h \
    src/solvers/COCG/COCG_MKL/COCG_Smooth_MKL.h \
    src/solvers/COCG/COCG_MKL/COCG_Di_Smooth_MKL.h \
    src/solvers/COCG/COCG_MKL/COCG_LLT_Smooth_MKL.h

DESTDIR = .
OBJECTS_DIR = build
MOC_DIR = build
RCC_DIR = build
UI_DIR = build
