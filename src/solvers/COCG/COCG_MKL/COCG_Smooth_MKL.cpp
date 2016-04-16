#if defined(_MSC_VER) && !defined(_CRT_SECURE_NO_WARNINGS)
#define _CRT_SECURE_NO_WARNINGS
#endif
#include "COCG_Smooth_MKL.h"
#include <cstdio>
#include <cmath>
#include <algorithm>

#if !defined(USE_OMP)
#include "../../../stubs/omp_stubs.h"
using namespace omp_stubs;
#else
#include <omp.h>
#endif

#if defined(_MSC_VER)
typedef long omp_int;
#else
typedef size_t omp_int;
#endif

// y := x
inline void normal_zcopy(size_t n, const complex<double> * x, complex<double> * y)
{
    MKL_INT m = (MKL_INT)n;
    MKL_Complex16 * x1 = reinterpret_cast<MKL_Complex16 *>(const_cast<complex<double> *>(x));
    MKL_Complex16 * y1 = reinterpret_cast<MKL_Complex16 *>(const_cast<complex<double> *>(y));
    cblas_zcopy(m, x1, 1, y1, 1);
}

// y := a * x + y
inline void normal_zaxpy(size_t n, const complex<double> & a, const complex<double> * x, complex<double> * y)
{
    MKL_INT m = (MKL_INT)n;
    MKL_Complex16 * a1 = reinterpret_cast<MKL_Complex16 *>(const_cast<complex<double> *>(&a));
    MKL_Complex16 * x1 = reinterpret_cast<MKL_Complex16 *>(const_cast<complex<double> *>(x));
    MKL_Complex16 * y1 = reinterpret_cast<MKL_Complex16 *>(const_cast<complex<double> *>(y));
    cblas_zaxpy(m, a1, x1, 1, y1, 1);
}

void COCG_Smooth_MKL::init(const size_t * gi_s, const size_t * gj_s, const complex<double> * di_s,
                          const complex<double> * gg_s, size_t n_s)
{
    gi = gi_s;
    gj = gj_s;
    di = di_s;
    gg = gg_s;
    n = n_s;

    delete [] r;
    delete [] z;
    delete [] p;
    delete [] s;
    delete [] xs;
    delete [] rs;

    r = new complex<double> [n];
    z = new complex<double> [n];
    p = new complex<double> [n];
    s = new complex<double> [n];
    xs = new complex<double> [n];
    rs = new complex<double> [n];

    delete [] ia;
    delete [] ja;
    delete [] aa;

    m = (MKL_INT)n;
    ia = new MKL_INT[n + 1];
    ia[0] = 1;
    for(MKL_INT i = 1; i <= m; i++)
        ia[i] = ia[i - 1] + (MKL_INT)(gi[i] - gi[i - 1]) + 1;
    ja = new MKL_INT[ia[n]];
    aa = new MKL_Complex16[ia[n]];
#pragma omp parallel for
    for(omp_int i = 0; i < (omp_int)n; i++)
    {
        size_t i_b = gi[i], i_e = gi[i + 1];
        for(size_t j = i_b; j < i_e; j++)
        {
            size_t ind = ia[i] + (j - gi[i]) - 1;
            ja[ind] = (MKL_INT)gj[j] + 1;
            aa[ind].real = gg[j].real();
            aa[ind].imag = gg[j].imag();
        }
        size_t ind = ia[i + 1] - 2;
        ja[ind] = i + 1;
        aa[ind].real = di[i].real();
        aa[ind].imag = di[i].imag();
    }
}

complex<double> COCG_Smooth_MKL::dot_prod_nocj(const complex<double> * a, const complex<double> * b) const
{
    complex<double> result;
    MKL_Complex16 * a1 = reinterpret_cast<MKL_Complex16 *>(const_cast<complex<double> *>(a));
    MKL_Complex16 * b1 = reinterpret_cast<MKL_Complex16 *>(const_cast<complex<double> *>(b));
    MKL_Complex16 * pres = reinterpret_cast<MKL_Complex16 *>(&result);
    cblas_zdotu_sub(m, a1, 1, b1, 1, pres);
    return result;
}

double COCG_Smooth_MKL::dot_prod_self(const complex<double> * a) const
{
    MKL_Complex16 * a1 = reinterpret_cast<MKL_Complex16 *>(const_cast<complex<double> *>(a));
    double result = cblas_dznrm2(m, a1, 1);
    return result * result;
}

double COCG_Smooth_MKL::dot_prod_real(const complex<double> * a, const complex<double> * b) const
{
    double d_p = 0.0;
#pragma omp parallel for reduction(+ : d_p)
    for(omp_int i = 0; i < (omp_int)n; i++)
        d_p += a[i].real() * b[i].real() + a[i].imag() * b[i].imag();
    return d_p;
}

void COCG_Smooth_MKL::mul_matrix(const complex<double> * f, complex<double> * x) const
{
    MKL_Complex16 * in_v = reinterpret_cast<MKL_Complex16 *>(const_cast<complex<double> *>(f));
    MKL_Complex16 * out_v = reinterpret_cast<MKL_Complex16 *>(const_cast<complex<double> *>(x));
    mkl_zcsrsymv("L", &m, aa, ia, ja, in_v, out_v);
}

bool COCG_Smooth_MKL::is_fpu_error(double x) const
{
    double y = x - x;
    return x != x || y != y;
}

void COCG_Smooth_MKL::solve(complex<double> * solution, const complex<double> * rp_s, double eps, size_t max_iter)
{
    eps *= eps;

    rp = rp_s;

    x0 = solution;
    mul_matrix(x0, r);

#pragma omp parallel for
    for(omp_int i = 0; i < (omp_int)n; i++)
    {
        r[i] = rp[i] - r[i];
        z[i] = r[i];
    }
    normal_zcopy(n, solution, xs);
    normal_zcopy(n, r, rs);
    normal_zcopy(n, z, p);

    complex<double> alpha, beta, prod_1, prod_2;
    double discr = 0.0, rp_norm, eta;

    rp_norm = dot_prod_self(rp);
    if(is_fpu_error(rp_norm))
    {
        fprintf(stderr, "Error: FPU error detected in right part!\n");
        return;
    }

    prod_1 = dot_prod_nocj(p, r);

    bool finished = false;

    size_t iter;
    for(iter = 0; iter <= max_iter && !finished; iter++)
    {
        discr = dot_prod_self(rs);
        if(is_fpu_error(discr))
        {
            fprintf(stderr, "Error: FPU error detected in (r, r)!\n");
            normal_zcopy(n, xs, solution);
            return;
        }

        double residual = discr / rp_norm;
        if(iter%10 == 0)
        {
            printf("COCG_Smooth_MKL [%d] Residual:\t%5lu\t%.3e\r", numThreads, (unsigned long)iter, sqrt(residual));
            fflush(stdout);
        }

        if(residual > eps)
        {
            mul_matrix(z, s);

            alpha = prod_1 / dot_prod_nocj(s, z);

            normal_zaxpy(n, alpha, z, x0);
            normal_zaxpy(n, -alpha, s, r);
#pragma omp parallel for
            for(omp_int i = 0; i < (omp_int)n ; i++)
            {
                p[i] = r[i];
                s[i] = r[i] - rs[i]; // r - s, сглаживатель
            }

            eta = - dot_prod_real(rs, s) / dot_prod_self(s);
            if(eta < 0.0) eta = 0.0;
            else if(eta > 1.0) eta = 1.0;
#pragma omp parallel for
            for(omp_int i = 0; i < (omp_int)n; i++)
            {
                double eta1 = 1.0 - eta;
                xs[i] = eta1 * xs[i] + eta * x0[i];
                rs[i] = eta1 * rs[i] + eta * r[i];
            }

            prod_2 = dot_prod_nocj(p, r);

            beta = prod_2 / prod_1;

            prod_1 = prod_2;

            normal_zaxpy(n, beta, z, p);
            swap(z, p);
        }
        else
            finished = true;
    }

//    mul_matrix(xs, r);
//#pragma omp parallel for
//            for(omp_int i = 0; i < (omp_int)n; i++)
//        r[i] = rp[i] - r[i];
//    discr = dot_prod_self(r);
    printf("COCG_Smooth_MKL [%d] Residual:\t%5lu\t%.3e\n", numThreads, (unsigned long)iter - 1, sqrt(discr / rp_norm));
    if(iter >= max_iter)
        printf("Soulution can`t found, iteration limit exceeded!\n");

    normal_zcopy(n, xs, solution);
}

COCG_Smooth_MKL::COCG_Smooth_MKL()
{
    numThreads = -1;
    if(numThreads <= 0)
    {
        char * env_num_threads = getenv("MKL_NUM_THREADS");
        if(env_num_threads)
            numThreads = atoi(env_num_threads);
    }
    if(numThreads <= 0)
    {
        char * env_num_threads = getenv("OMP_NUM_THREADS");
        if(env_num_threads)
            numThreads = atoi(env_num_threads);
    }
    if(numThreads <= 0)
    {
        numThreads = 1;
    }
    mkl_set_num_threads(numThreads);
    omp_set_num_threads(numThreads);
    numThreads = std::max(mkl_get_max_threads(), omp_get_max_threads());
    r = x0 = z = p = s = xs = rs = NULL;
    ia = ja = NULL;
    aa = NULL;
}

COCG_Smooth_MKL::~COCG_Smooth_MKL()
{
    delete [] r;
    delete [] z;
    delete [] p;
    delete [] s;
    delete [] xs;
    delete [] rs;
    delete [] ia;
    delete [] ja;
    delete [] aa;
}

// =================================================================================================

#if defined(USE_MKL) && defined(_MSC_VER)

//#define MKL_TBB_THREADS
//#define MKL_OMP_THREADS

#if defined(_WIN32)
#if !defined(_WIN64)
#pragma comment(lib, "mkl_intel_c.lib")
#else
#pragma comment(lib, "mkl_intel_lp64.lib")
#endif
#pragma comment(lib, "mkl_core.lib")
#if defined(MKL_TBB_THREADS)
#pragma comment(lib, "mkl_tbb_thread.lib")
#pragma comment(lib, "tbb.lib")
#elif defined(MKL_OMP_THREADS)
#pragma comment(lib, "mkl_intel_thread.lib")
#pragma comment(lib, "libiomp5md.lib")
#else
#pragma comment(lib, "mkl_sequential.lib")
#endif
#endif

#endif
