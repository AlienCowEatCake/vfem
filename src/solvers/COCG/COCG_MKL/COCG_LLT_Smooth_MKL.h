#if !defined COCG_LLT_SMOOTH_MKL_H_INCLUDED
#define COCG_LLT_SMOOTH_MKL_H_INCLUDED

//#define USE_MKL

#include <cstdlib>
#include <complex>
#include <omp.h>

#if !defined USE_MKL
#if !defined MKL_INT
#define MKL_INT int
typedef struct _MKL_Complex16 { double real; double imag; } MKL_Complex16;
void cblas_zcopy(const MKL_INT N, const void * X, const MKL_INT incX, void * Y, const MKL_INT incY);
void cblas_zaxpy(const MKL_INT N, const void * alpha, const void * X, const MKL_INT incX, void * Y, const MKL_INT incY);
double cblas_dznrm2(const MKL_INT N, const void * X, const MKL_INT incX);
void cblas_zdotu_sub(const MKL_INT N, const void * X, const MKL_INT incX, const void * Y, const MKL_INT incY, void * dotc);
void mkl_zcsrsymv(const char * uplo, const MKL_INT * m, const MKL_Complex16 * a, const MKL_INT * ia,  const MKL_INT * ja, const MKL_Complex16 * x,  MKL_Complex16 * y);
void mkl_zcsrtrsv(const char * uplo, const char * transa, const char * diag, const MKL_INT * m, const MKL_Complex16 * a, const MKL_INT * ia, const MKL_INT * ja, const MKL_Complex16 * x, MKL_Complex16 * y);
void mkl_set_num_threads(int n);
#endif
#else
#include <mkl.h>
#include <mkl_blas.h>
#include <mkl_cblas.h>
#include <mkl_spblas.h>
#endif

#include "../../solver_interface.h"

using namespace std;

class COCG_LLT_Smooth_MKL : public solver_interface<complex<double>, size_t>
{
public:
    void init(const size_t * gi_s, const size_t * gj_s, const complex<double> * di_s,
              const complex<double> * gg_s, size_t n_s);
    void solve(complex<double> * solution, const complex<double> * rp_s, double eps, size_t max_iter);

    COCG_LLT_Smooth_MKL();
    ~COCG_LLT_Smooth_MKL();
private:
    void make_LLT_decomposition();
    void mul_matrix(const complex<double> * f, complex<double> * x) const;
    void solve_LLT(const complex<double> * f, complex<double> * x) const;
    complex<double> dot_prod_nocj(const complex<double> * a, const complex<double> * b) const;
    double dot_prod_self(const complex<double> * a) const;
    double dot_prod_real(const complex<double> * a, const complex<double> * b) const;
    bool is_fpu_error(double x) const;

    size_t n;
    const size_t * gi, * gj;
    const complex<double> * di, * gg, * rp;
    complex<double> * r, * x0, * z, * p, * s, * xs, * rs;
    MKL_Complex16 * L_aa;
    mutable MKL_Complex16 * LLT_tmp;

    int numThreads;
    MKL_INT m;
    MKL_INT * ia, * ja;
    MKL_Complex16 * aa;
};

#endif // COCG_LLT_SMOOTH_MKL_H_INCLUDED
