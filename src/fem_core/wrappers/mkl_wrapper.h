#if !defined(WRAPPERS_MKL_WRAPPER_H_INCLUDED)
#define WRAPPERS_MKL_WRAPPER_H_INCLUDED

// *************************************************************************************************

#if defined(USE_MKL)

#include <mkl.h>
#include <mkl_blas.h>
#include <mkl_cblas.h>
#include <mkl_spblas.h>

#else

namespace fem_core { namespace wrappers { namespace mkl_stubs {

typedef int MKL_INT;
typedef struct _MKL_Complex8 { float real; float imag; } MKL_Complex8;
typedef struct _MKL_Complex16 { double real; double imag; } MKL_Complex16;

void cblas_scopy(const MKL_INT N, const void * X, const MKL_INT incX, void * Y, const MKL_INT incY);
void cblas_dcopy(const MKL_INT N, const void * X, const MKL_INT incX, void * Y, const MKL_INT incY);
void cblas_ccopy(const MKL_INT N, const void * X, const MKL_INT incX, void * Y, const MKL_INT incY);
void cblas_zcopy(const MKL_INT N, const void * X, const MKL_INT incX, void * Y, const MKL_INT incY);

void cblas_zaxpy(const MKL_INT N, const void * alpha, const void * X, const MKL_INT incX, void * Y, const MKL_INT incY);

double cblas_dznrm2(const MKL_INT N, const void * X, const MKL_INT incX);

void cblas_zdotu_sub(const MKL_INT N, const void * X, const MKL_INT incX, const void * Y, const MKL_INT incY, void * dotc);

void cblas_zdotc_sub(const MKL_INT N, const void * X, const MKL_INT incX, const void * Y, const MKL_INT incY, void * dotc);

void mkl_zcsrsymv(const char * uplo, const MKL_INT * m, const MKL_Complex16 * a, const MKL_INT * ia,  const MKL_INT * ja, const MKL_Complex16 * x,  MKL_Complex16 * y);

void mkl_zcsrtrsv(const char * uplo, const char * transa, const char * diag, const MKL_INT * m, const MKL_Complex16 * a, const MKL_INT * ia, const MKL_INT * ja, const MKL_Complex16 * x, MKL_Complex16 * y);

void mkl_set_num_threads(int n);

int mkl_get_max_threads();

}}} // namespace fem_core::wrappers::mkl_stubs

using namespace fem_core::wrappers::mkl_stubs;

#endif

// *************************************************************************************************

namespace fem_core { namespace wrappers { namespace mkl {

/**
 * @brief Установить количество тредов для MKL равным переменной окружения MKL_NUM_THREADS
 */
void wrapper_mkl_set_env_max_threads();

}}} // namespace fem_core::wrappers::mkl

// *************************************************************************************************

#endif // WRAPPERS_MKL_WRAPPER_H_INCLUDED
