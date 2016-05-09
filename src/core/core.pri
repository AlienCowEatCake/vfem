INCLUDEPATH += $$PWD
DEPENDPATH += $$PWD
CONFIG += object_with_source object_parallel_to_source
CONFIG += no_batch

HEADERS += \
    $$PWD/containers/generic/array_t.h \
    $$PWD/containers/generic/matrix_t.h \
    $$PWD/containers/generic/trio.h \
    $$PWD/containers/geometry/edge_basic.h \
    $$PWD/containers/geometry/face_triangle_basic.h \
    $$PWD/containers/geometry/point3_t.h \
    $$PWD/containers/geometry/vector3_t.h \
    $$PWD/containers/tree/octree.h \
    $$PWD/containers/tree/quadtree.h \
    $$PWD/evaluator/evaluator_internal/jit/common.h \
    $$PWD/evaluator/evaluator_internal/jit/compile_extcall.h \
    $$PWD/evaluator/evaluator_internal/jit/compile_inline.h \
    $$PWD/evaluator/evaluator_internal/jit/complex_templates.h \
    $$PWD/evaluator/evaluator_internal/jit/func_templates.h \
    $$PWD/evaluator/evaluator_internal/jit/opcodes.h \
    $$PWD/evaluator/evaluator_internal/jit/oper_templates.h \
    $$PWD/evaluator/evaluator_internal/jit/real_templates.h \
    $$PWD/evaluator/evaluator_internal/calculate.h \
    $$PWD/evaluator/evaluator_internal/evaluator_object.h \
    $$PWD/evaluator/evaluator_internal/misc.h \
    $$PWD/evaluator/evaluator_internal/parse.h \
    $$PWD/evaluator/evaluator_internal/simplify.h \
    $$PWD/evaluator/evaluator_internal/transition_table.h \
    $$PWD/evaluator/evaluator_internal/type_detection.h \
    $$PWD/evaluator/evaluator_internal/var_container.h \
    $$PWD/evaluator/evaluator.h \
    $$PWD/evaluator/evaluator_operations.h \
    $$PWD/evaluator/evaluator_xyz.h \
    $$PWD/solvers/CSRC/preconditioners/preconditioner_interface.h \
    $$PWD/solvers/CSRC/preconditioners/Nothing/preconditioner_Nothing.h \
    $$PWD/solvers/CSRC/preconditioners/Nothing/preconditioner_Nothing_MKL.h \
    $$PWD/solvers/CSRC/preconditioners/Nothing/preconditioner_Nothing_OpenMP.h \
    $$PWD/solvers/CSRC/preconditioners/Di/preconditioner_Di.h \
    $$PWD/solvers/CSRC/preconditioners/Di/preconditioner_Di_MKL.h \
    $$PWD/solvers/CSRC/preconditioners/Di/preconditioner_Di_OpenMP.h \
    $$PWD/solvers/CSRC/preconditioners/GS/preconditioner_GS.h \
    $$PWD/solvers/CSRC/preconditioners/LDLT/preconditioner_LDLT.h \
    $$PWD/solvers/CSRC/preconditioners/LLT/preconditioner_LLT.h \
    $$PWD/solvers/CSRC/preconditioners/LLT/preconditioner_LLT_MKL.h \
    $$PWD/solvers/CSRC/symmetric/symmetric_solver_interface.h \
    $$PWD/solvers/CSRC/symmetric/complex/COCG/COCG.h \
    $$PWD/solvers/CSRC/symmetric/complex/COCG/COCG_Smooth.h \
    $$PWD/solvers/CSRC/symmetric/complex/COCG/COCG_Smooth_MKL.h \
    $$PWD/solvers/CSRC/symmetric/complex/COCG/COCG_Smooth_OpenMP.h \
    $$PWD/solvers/CSRC/symmetric/complex/COCR/COCR.h \
    $$PWD/solvers/CSRC/symmetric/complex/COCR/COCR_Smooth.h \
    $$PWD/solvers/CSRC/symmetric/complex/BiCG_Complex/BiCG_Complex.h \
    $$PWD/solvers/CSRC/symmetric/complex/BiCG_Complex/BiCG_Complex_Smooth.h \
    $$PWD/solvers/CSRC/symmetric/complex/BiCGStab_Complex/BiCGStab_Complex.h \
    $$PWD/solvers/CSRC/symmetric/complex/BiCGStab_Complex/BiCGStab_Complex_Smooth.h \
    $$PWD/solvers/CSRC/symmetric/complex/GMRES_Complex/GMRES_Complex.h \
    $$PWD/solvers/CSRC/symmetric/complex/GMRES_Complex/GMRES_Complex_MKL.h \
    $$PWD/solvers/CSRC/symmetric/complex/GMRES_Complex/GMRES_Complex_OpenMP.h \
    $$PWD/solvers/CSRC/solvers_factory.h \
    $$PWD/utils/cxxversion.h \
    $$PWD/utils/fpu.h \
    $$PWD/utils/inifile.h \
    $$PWD/utils/progress.h \
    $$PWD/utils/strings.h \
    $$PWD/utils/timers.h \
    $$PWD/wrappers/mkl_wrapper.h \
    $$PWD/wrappers/omp_wrapper.h

SOURCES += \
    $$PWD/evaluator/evaluator_internal/jit/common.cpp \
    $$PWD/evaluator/evaluator_internal/jit/func_templates.cpp \
    $$PWD/evaluator/evaluator_internal/jit/oper_templates.cpp \
    $$PWD/evaluator/evaluator_internal/jit/real_templates.cpp \
    $$PWD/evaluator/evaluator_internal/transition_table.cpp \
    $$PWD/solvers/CSRC/symmetric/complex/COCG/COCG.cpp \
    $$PWD/solvers/CSRC/symmetric/complex/COCG/COCG_Smooth.cpp \
    $$PWD/solvers/CSRC/symmetric/complex/COCG/COCG_Smooth_MKL.cpp \
    $$PWD/solvers/CSRC/symmetric/complex/COCG/COCG_Smooth_OpenMP.cpp \
    $$PWD/solvers/CSRC/symmetric/complex/COCR/COCR.cpp \
    $$PWD/solvers/CSRC/symmetric/complex/COCR/COCR_Smooth.cpp \
    $$PWD/solvers/CSRC/symmetric/complex/BiCG_Complex/BiCG_Complex.cpp \
    $$PWD/solvers/CSRC/symmetric/complex/BiCG_Complex/BiCG_Complex_Smooth.cpp \
    $$PWD/solvers/CSRC/symmetric/complex/BiCGStab_Complex/BiCGStab_Complex.cpp \
    $$PWD/solvers/CSRC/symmetric/complex/BiCGStab_Complex/BiCGStab_Complex_Smooth.cpp \
    $$PWD/solvers/CSRC/symmetric/complex/GMRES_Complex/GMRES_Complex.cpp \
    $$PWD/solvers/CSRC/symmetric/complex/GMRES_Complex/GMRES_Complex_MKL.cpp \
    $$PWD/solvers/CSRC/symmetric/complex/GMRES_Complex/GMRES_Complex_OpenMP.cpp \
    $$PWD/solvers/CSRC/solvers_factory.cpp \
    $$PWD/utils/inifile.cpp \
    $$PWD/utils/progress.cpp \
    $$PWD/utils/strings.cpp \
    $$PWD/utils/timers.cpp \
    $$PWD/wrappers/mkl_wrapper.cpp \
    $$PWD/wrappers/omp_wrapper.cpp

unix:!macx:QMAKE_LIBS += -lrt

*g++*|*clang* {
    QMAKE_CXXFLAGS_RELEASE -= -O2
    QMAKE_CXXFLAGS_RELEASE *= -O3
    QMAKE_CXXFLAGS_RELEASE *= -DNDEBUG
}

*msvc* {
    QMAKE_CXXFLAGS_RELEASE -= -O2
    QMAKE_CXXFLAGS_RELEASE *= -Ox
    QMAKE_CXXFLAGS_RELEASE -= -GS
    QMAKE_CXXFLAGS_RELEASE *= -GS-
}

use_mkl {
    win32-msvc201*|win32-icc* {
        contains(QMAKE_HOST.arch, x86_64) {
            PROGRAMFILES = $$system("echo %ProgramFiles(x86)%")
        } else {
            PROGRAMFILES = $$(ProgramFiles)
        }
        MKLROOT =
        KNOWN_MKLROOTS = \
            $${PROGRAMFILES}/IntelSWTools/compilers_and_libraries_2016.0.110/windows/mkl \
            $${PROGRAMFILES}/IntelSWTools/compilers_and_libraries_2016.1.146/windows/mkl \
            $${PROGRAMFILES}/IntelSWTools/compilers_and_libraries_2016.2.180/windows/mkl
        for(TEST_MKLROOT, KNOWN_MKLROOTS):exists($${TEST_MKLROOT}):MKLROOT=$${TEST_MKLROOT}
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
        MKLROOT =
        KNOWN_MKLROOTS = \
            /opt/intel/compilers_and_libraries_2016.0.109/linux/mkl \
            /opt/intel/compilers_and_libraries_2016.1.150/linux/mkl \
            /opt/intel/compilers_and_libraries_2016.2.181/linux/mkl
        for(TEST_MKLROOT, KNOWN_MKLROOTS):exists($${TEST_MKLROOT}):MKLROOT=$${TEST_MKLROOT}
        exists($${MKLROOT}) {
            MKLLIB = $${MKLROOT}/lib/intel64
            QMAKE_CXXFLAGS  -= -ansi
            QMAKE_CXXFLAGS  *= -I$${MKLROOT}/include -DUSE_MKL -std=c++0x
            QMAKE_LFLAGS    *= -Wl,-rpath -Wl,$${MKLLIB} -Wl,--no-as-needed -L$${MKLLIB}
            QMAKE_LIBS      *= -lmkl_intel_lp64 -lmkl_core -lmkl_gnu_thread
            #QMAKE_LIBS      *= -Wl,--start-group $${MKLLIB}/libmkl_intel_lp64.a $${MKLLIB}/libmkl_core.a $${MKLLIB}/libmkl_gnu_thread.a -Wl,--end-group
            QMAKE_LIBS      *= -ldl -lpthread -lm
            !use_omp:CONFIG+=use_omp
        }
    }
}

use_omp {
    *g++* {
        QMAKE_CXXFLAGS += -fopenmp
        QMAKE_LFLAGS += -fopenmp
        DEFINES += USE_OMP
    }
    *clang* {
        # Clang пока что не умеет OpenMP
        QMAKE_CXXFLAGS += -Wno-unknown-pragmas
    }
    *msvc* {
        QMAKE_CXXFLAGS += -openmp
        DEFINES += USE_OMP
    }
} else {
    *g++*|*clang* {
        QMAKE_CXXFLAGS += -Wno-unknown-pragmas
    }
}
