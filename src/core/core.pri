INCLUDEPATH += $$PWD
DEPENDPATH += $$PWD
CONFIG += object_with_source object_parallel_to_source
CONFIG += no_batch

HEADERS += \
    $$files($$PWD/containers/fem/*.h) \
    $$files($$PWD/containers/generic/*.h) \
    $$files($$PWD/containers/geometry/*.h) \
    $$files($$PWD/containers/tree/*.h) \
    $$files($$PWD/cubatures/*.h) \
    $$files($$PWD/evaluator/evaluator_internal/jit/*.h) \
    $$files($$PWD/evaluator/evaluator_internal/*.h) \
    $$files($$PWD/evaluator/*.h) \
    $$files($$PWD/solvers/CSLR/preconditioners/Nothing/*.h) \
    $$files($$PWD/solvers/CSLR/preconditioners/Di/*.h) \
    $$files($$PWD/solvers/CSLR/preconditioners/GS/*.h) \
    $$files($$PWD/solvers/CSLR/preconditioners/LDLT/*.h) \
    $$files($$PWD/solvers/CSLR/preconditioners/LLT/*.h) \
    $$files($$PWD/solvers/CSLR/preconditioners/*.h) \
    $$files($$PWD/solvers/CSLR/symmetric/complex/COCG/*.h) \
    $$files($$PWD/solvers/CSLR/symmetric/complex/COCR/*.h) \
    $$files($$PWD/solvers/CSLR/symmetric/complex/BiCG_Complex/*.h) \
    $$files($$PWD/solvers/CSLR/symmetric/complex/BiCGStab_Complex/*.h) \
    $$files($$PWD/solvers/CSLR/symmetric/complex/GMRES_Complex/*.h) \
    $$files($$PWD/solvers/CSLR/symmetric/*.h) \
    $$files($$PWD/solvers/CSLR/*.h) \
    $$files($$PWD/utils/*.h) \
    $$files($$PWD/wrappers/*.h)

SOURCES += \
    $$files($$PWD/cubatures/*.cpp) \
    $$files($$PWD/evaluator/evaluator_internal/jit/*.cpp) \
    $$files($$PWD/evaluator/evaluator_internal/*.cpp) \
    $$files($$PWD/solvers/CSLR/symmetric/complex/COCG/*.cpp) \
    $$files($$PWD/solvers/CSLR/symmetric/complex/COCR/*.cpp) \
    $$files($$PWD/solvers/CSLR/symmetric/complex/BiCG_Complex/*.cpp) \
    $$files($$PWD/solvers/CSLR/symmetric/complex/BiCGStab_Complex/*.cpp) \
    $$files($$PWD/solvers/CSLR/symmetric/complex/GMRES_Complex/*.cpp) \
    $$files($$PWD/solvers/CSLR/*.cpp) \
    $$files($$PWD/utils/*.cpp) \
    $$files($$PWD/wrappers/*.cpp)

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
            $${PROGRAMFILES}/IntelSWTools/compilers_and_libraries_2016.2.180/windows/mkl \
            $${PROGRAMFILES}/IntelSWTools/compilers_and_libraries_2016.3.207/windows/mkl
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
    linux-g++* {
        MKLROOT =
        KNOWN_MKLROOTS = \
            /opt/intel/compilers_and_libraries_2016.0.109/linux/mkl \
            /opt/intel/compilers_and_libraries_2016.1.150/linux/mkl \
            /opt/intel/compilers_and_libraries_2016.2.181/linux/mkl \
            /opt/intel/compilers_and_libraries_2016.3.210/linux/mkl
        for(TEST_MKLROOT, KNOWN_MKLROOTS):exists($${TEST_MKLROOT}):MKLROOT=$${TEST_MKLROOT}
        exists($${MKLROOT}) {
            MKLLIBDIR =
            MKLLIBS =
            linux-g++*-64 {
                MKLLIBDIR = $${MKLROOT}/lib/intel64
                MKLLIBS = -lmkl_intel_lp64 -lmkl_core -lmkl_gnu_thread
                #MKLLIBS = -Wl,--start-group $${MKLLIBDIR}/libmkl_intel_lp64.a $${MKLLIBDIR}/libmkl_core.a $${MKLLIBDIR}/libmkl_gnu_thread.a -Wl,--end-group
            }
            linux-g++*-32 {
                MKLLIBDIR = $${MKLROOT}/lib/ia32
                MKLLIBS = -lmkl_intel -lmkl_core -lmkl_gnu_thread
                #MKLLIBS = -Wl,--start-group $${MKLLIBDIR}/libmkl_intel.a $${MKLLIBDIR}/libmkl_core.a $${MKLLIBDIR}/libmkl_gnu_thread.a -Wl,--end-group
            }
            exists($${MKLLIBDIR}) {
                QMAKE_CXXFLAGS  *= -I$${MKLROOT}/include -DUSE_MKL -Wno-long-long
                QMAKE_LFLAGS    *= -Wl,-rpath -Wl,$${MKLLIBDIR} -Wl,--no-as-needed -L$${MKLLIBDIR}
                QMAKE_LIBS      *= $${MKLLIBS} -ldl -lpthread -lm
                !use_omp:CONFIG+=use_omp
            }
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
