#ifndef MKL_STUBS_H
#define MKL_STUBS_H

#if !defined(USE_OMP)
#include "omp_stubs.h"
using namespace omp_stubs;
#else
#include <omp.h>
#endif

namespace mkl_stubs
{

typedef int MKL_INT;
typedef struct _MKL_Complex16 { double real; double imag; } MKL_Complex16;

void cblas_zcopy(const MKL_INT N, const void * X, const MKL_INT incX, void * Y, const MKL_INT incY);

void cblas_zaxpy(const MKL_INT N, const void * alpha, const void * X, const MKL_INT incX, void * Y, const MKL_INT incY);

double cblas_dznrm2(const MKL_INT N, const void * X, const MKL_INT incX);

void cblas_zdotu_sub(const MKL_INT N, const void * X, const MKL_INT incX, const void * Y, const MKL_INT incY, void * dotc);

void mkl_zcsrsymv(const char * uplo, const MKL_INT * m, const MKL_Complex16 * a, const MKL_INT * ia,  const MKL_INT * ja, const MKL_Complex16 * x,  MKL_Complex16 * y);

void mkl_zcsrtrsv(const char * uplo, const char * transa, const char * diag, const MKL_INT * m, const MKL_Complex16 * a, const MKL_INT * ia, const MKL_INT * ja, const MKL_Complex16 * x, MKL_Complex16 * y);

void mkl_set_num_threads(int n);

int mkl_get_max_threads();

}

#endif // MKL_STUBS_H
