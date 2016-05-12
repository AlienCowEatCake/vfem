INCLUDEPATH += $$PWD
DEPENDPATH += $$PWD
CONFIG += object_with_source object_parallel_to_source
CONFIG += no_batch

HEADERS += \
    $$PWD/containers/generic/*.h \
    $$PWD/containers/geometry/*.h \
    $$PWD/containers/tree/*.h \
    $$PWD/evaluator/evaluator_internal/jit/*.h \
    $$PWD/evaluator/evaluator_internal/*.h \
    $$PWD/evaluator/*.h \
    $$PWD/solvers/CSRC/preconditioners/Nothing/*.h \
    $$PWD/solvers/CSRC/preconditioners/Di/*.h \
    $$PWD/solvers/CSRC/preconditioners/GS/*.h \
    $$PWD/solvers/CSRC/preconditioners/LDLT/*.h \
    $$PWD/solvers/CSRC/preconditioners/LLT/*.h \
    $$PWD/solvers/CSRC/preconditioners/*.h \
    $$PWD/solvers/CSRC/symmetric/complex/COCG/*.h \
    $$PWD/solvers/CSRC/symmetric/complex/COCR/*.h \
    $$PWD/solvers/CSRC/symmetric/complex/BiCG_Complex/*.h \
    $$PWD/solvers/CSRC/symmetric/complex/BiCGStab_Complex/*.h \
    $$PWD/solvers/CSRC/symmetric/complex/GMRES_Complex/*.h \
    $$PWD/solvers/CSRC/symmetric/*.h \
    $$PWD/solvers/CSRC/*.h \
    $$PWD/utils/*.h \
    $$PWD/wrappers/*.h

SOURCES += \
    $$PWD/evaluator/evaluator_internal/jit/*.cpp \
    $$PWD/evaluator/evaluator_internal/*.cpp \
    $$PWD/solvers/CSRC/symmetric/complex/COCG/*.cpp \
    $$PWD/solvers/CSRC/symmetric/complex/COCR/*.cpp \
    $$PWD/solvers/CSRC/symmetric/complex/BiCG_Complex/*.cpp \
    $$PWD/solvers/CSRC/symmetric/complex/BiCGStab_Complex/*.cpp \
    $$PWD/solvers/CSRC/symmetric/complex/GMRES_Complex/*.cpp \
    $$PWD/solvers/CSRC/*.cpp \
    $$PWD/utils/*.cpp \
    $$PWD/wrappers/*.cpp

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
            QMAKE_CXXFLAGS  *= -I$${MKLROOT}/include -DUSE_MKL -Wno-long-long
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
