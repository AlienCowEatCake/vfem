#if !defined(GMRES_COMPLEX_LLT_MKL_H_INCLUDED)
#define GMRES_COMPLEX_LLT_MKL_H_INCLUDED

#include <cstdlib>
#include <complex>

#if !defined(USE_MKL)
#include "../../../stubs/mkl_stubs.h"
using namespace mkl_stubs;
#else
#include <mkl.h>
#include <mkl_blas.h>
#include <mkl_cblas.h>
#include <mkl_spblas.h>
#endif

#include "../../solver_interface.h"

using namespace std;

class GMRES_Complex_LLT_MKL : public solver_interface<complex<double>, size_t>
{
public:
    void init(const size_t * gi_s, const size_t * gj_s, const complex<double> * di_s,
              const complex<double> * gg_s, size_t n_s);
    void solve(complex<double> * solution, const complex<double> * rp_s,
               double eps, size_t max_iter);

    GMRES_Complex_LLT_MKL();
    ~GMRES_Complex_LLT_MKL();
private:
    void make_LLT_decomposition();
    void mul_matrix(const complex<double> * f, complex<double> * x) const;
    void solve_L(const complex<double> * f, complex<double> * x) const;
    void solve_LT(complex<double> * f, complex<double> * x) const;
    void mul_LT(const complex<double> * f, complex<double> * x) const;
    void solve_LTAL(const complex<double> * f, complex<double> * x, complex<double> * tmp) const;
    complex<double> dot_prod(const complex<double> * a, const complex<double> * b) const;
    double dot_prod_self(const complex<double> * a) const;
    bool is_fpu_error(double x) const;
    void copy_vec(const complex<double> * f, complex<double> * x) const;

    size_t n;
    const size_t * gi, * gj;
    const complex<double> * di, * gg, * rp;
    complex<double> * r, * x0;

    size_t m, m_curr;
    complex<double> * t, * w, ** VT, ** H, * d;

    MKL_Complex16 * L_aa;

    int numThreads;
    MKL_INT m_mkl;
    MKL_INT * ia, * ja;
    MKL_Complex16 * aa;
};

#endif // GMRES_COMPLEX_LLT_MKL_H_INCLUDED
