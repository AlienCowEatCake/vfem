#if defined(_MSC_VER) && !defined(_CRT_SECURE_NO_WARNINGS)
#define _CRT_SECURE_NO_WARNINGS
#endif
#include "GMRES_Complex_MKL.h"
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

// y := a * x + y
inline void normal_zaxpy(size_t n, const complex<double> & a, const complex<double> * x, complex<double> * y)
{
    MKL_INT m = (MKL_INT)n;
    MKL_Complex16 * a1 = reinterpret_cast<MKL_Complex16 *>(const_cast<complex<double> *>(&a));
    MKL_Complex16 * x1 = reinterpret_cast<MKL_Complex16 *>(const_cast<complex<double> *>(x));
    MKL_Complex16 * y1 = reinterpret_cast<MKL_Complex16 *>(const_cast<complex<double> *>(y));
    cblas_zaxpy(m, a1, x1, 1, y1, 1);
}

void GMRES_Complex_MKL::init(const size_t * gi_s, const size_t * gj_s, const complex<double> * di_s,
                             const complex<double> * gg_s, size_t n_s)
{
    gi = gi_s;
    gj = gj_s;
    di = di_s;
    gg = gg_s;
    n = n_s;

    delete [] r;
    delete [] w;
    delete [] d;
    if(H && m)
        for(size_t i = 0; i <= m; i++)
            delete [] H[i];
    delete [] H;
    if(VT && m)
        for(size_t i = 0; i < m; i++)
            delete [] VT[i];
    delete [] VT;

    r = new complex<double> [n];
    w = new complex<double> [n];

    // Глубина метода
    m = m_curr = 5;
//    FILE * fp = fopen("gmres_m.txt", "r");
//    if(fp)
//    {
//        unsigned int _m;
//        fscanf(fp, "%u", & _m);
//        m = m_curr = _m;
//        fclose(fp);
//    }
    char * env_m = getenv("GMRES_M");
    if(env_m)
    {
        int temp_m = atoi(env_m);
        if(temp_m >= 1)
            m = m_curr = (size_t)temp_m;
    }

    VT = new complex<double> * [m];
    for(size_t i = 0; i < m; i++)
        VT[i] = new complex<double> [n];
    H = new complex<double> * [m + 1];
    for(size_t i = 0; i < m + 1; i++)
    {
        H[i] = new complex<double> [m];
        for(size_t j = 0; j < m; j++)
            H[i][j] = 0;
    }
    d = new complex<double> [m + 1];

    m_mkl = (MKL_INT)n;
    ia = new MKL_INT[n + 1];
    ia[0] = 1;
    for(MKL_INT i = 1; i <= m_mkl; i++)
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

complex<double> GMRES_Complex_MKL::dot_prod(const complex<double> * a, const complex<double> * b) const
{
    complex<double> result;
    MKL_Complex16 * a1 = reinterpret_cast<MKL_Complex16 *>(const_cast<complex<double> *>(a));
    MKL_Complex16 * b1 = reinterpret_cast<MKL_Complex16 *>(const_cast<complex<double> *>(b));
    MKL_Complex16 * pres = reinterpret_cast<MKL_Complex16 *>(&result);
    cblas_zdotc_sub(m_mkl, a1, 1, b1, 1, pres);
    return result;
}

double GMRES_Complex_MKL::dot_prod_self(const complex<double> * a) const
{
    MKL_Complex16 * a1 = reinterpret_cast<MKL_Complex16 *>(const_cast<complex<double> *>(a));
    double result = cblas_dznrm2(m_mkl, a1, 1);
    return result * result;
}

void GMRES_Complex_MKL::mul_matrix(const complex<double> * f, complex<double> * x) const
{
    MKL_Complex16 * in_v = reinterpret_cast<MKL_Complex16 *>(const_cast<complex<double> *>(f));
    MKL_Complex16 * out_v = reinterpret_cast<MKL_Complex16 *>(const_cast<complex<double> *>(x));
    mkl_zcsrsymv("L", &m_mkl, aa, ia, ja, in_v, out_v);
}

void GMRES_Complex_MKL::copy_vec(const complex<double> * f, complex<double> * x) const
{
    MKL_Complex16 * x1 = reinterpret_cast<MKL_Complex16 *>(const_cast<complex<double> *>(f));
    MKL_Complex16 * y1 = reinterpret_cast<MKL_Complex16 *>(const_cast<complex<double> *>(x));
    cblas_zcopy(m_mkl, x1, 1, y1, 1);
}

bool GMRES_Complex_MKL::is_fpu_error(double x) const
{
    double y = x - x;
    return x != x || y != y;
}

void GMRES_Complex_MKL::solve(complex<double> * solution, const complex<double> * rp_s,
                              double eps, size_t max_iter)
{
    rp = rp_s;

    x0 = solution;

    mul_matrix(solution, r);
#pragma omp parallel for
    for(omp_int i = 0; i < (omp_int)n ; i++)
        r[i] = rp[i] - r[i];

    double discr = 0.0, rp_norm = sqrt(dot_prod_self(rp));
    double residual = 0.0, residual_prev;

    if(is_fpu_error(rp_norm))
    {
        fprintf(stderr, "Error: FPU error detected in right part!\n");
        return;
    }

    bool finished = false;

    size_t iter;
    for(iter = 0; iter <= max_iter && !finished; iter++)
    {
        discr = sqrt(dot_prod_self(r));
        if(is_fpu_error(discr))
        {
            fprintf(stderr, "Error: FPU error detected in (r, r)!\n");
            return;
        }

        residual_prev = residual;
        residual = discr / rp_norm;
        //if(iter%10 == 0)
        {
            printf("GMRES_Complex_MKL(%lu) [%d] Residual:\t%5lu\t%.3e\r", (unsigned long)m, numThreads, (unsigned long)iter, residual);
            fflush(stdout);
        }

        if(residual > eps && fabs(residual - residual_prev) / (residual) > 1e-15)
        {
            // V
            #pragma omp parallel for
            for(omp_int i = 0; i < (omp_int)n; i++)
                VT[0][i] = r[i] / discr;

            // d
            d[0] = discr;
            for(size_t i = 1; i <= m_curr; i++)
                d[i] = 0.0;

            for(size_t mu = 0; mu < m_curr; mu++)
            {
                // w
                mul_matrix(VT[mu], w);

                // H_lambda_mu
                for(size_t lambda = 0; lambda <= mu; lambda++)
                    H[lambda][mu] = dot_prod(VT[lambda], w);

                // tilde_v_mu+1
                for(size_t lambda = 0; lambda <= mu; lambda++)
                {
                    #pragma omp parallel for
                    for(omp_int i = 0; i < (omp_int)n; i++)
                        w[i] -= VT[lambda][i] * H[lambda][mu];
                }

                // H_mu+1_mu
                double vnorm = sqrt(dot_prod_self(w));
                H[mu + 1][mu] = vnorm;

                if(abs(vnorm) < 1e-15)
                {
                    if(mu > 1)  m_curr = mu;
                    else        m_curr = 1;
                    break;
                }

                // v_mu+1
                if(mu + 1 < m_curr)
                {
                    #pragma omp parallel for
                    for(omp_int i = 0; i < (omp_int)n; i++)
                        VT[mu + 1][i] = w[i] / vnorm;
                }
            }

            // H(i), d(i)
            for(size_t i = 0; i < m_curr; i++)
            {
                double t1 = H[i + 1][i].real();
                double t2 = abs(H[i][i]);
                double t3 = 1.0 / sqrt(t2 * t2 + t1 * t1);
                double s = t1 * t3;
                complex<double> c = H[i][i] * t3;

                for(size_t j = i; j < m_curr; j++)
                {
                    complex<double> H_ij = H[i][j];
                    H[i][j] = H_ij * conj(c) + H[i + 1][j] * s;
                    H[i + 1][j] = - H_ij * s + H[i + 1][j] * c;
                }
                complex<double> d_i = d[i];
                d[i] = d_i * conj(c) + d[i + 1] * s;
                d[i + 1] = - d_i * s + d[i + 1] * c;
            }

            // Hz = d
            for(size_t i = 0; i < m_curr; i++)
            {
                d[i] = d[i] / H[i][i];
                for(size_t j = i + 1; j < m_curr; j++)
                    H[i][j] /= H[i][i];
            }
            for(size_t i = m_curr - 1; i > 0; i--)
                for(size_t j = i - 1, ju = i; ju > 0; j--, ju--)
                    d[j] -= H[j][i] * d[i];

            // x = x + Vz
            #pragma omp parallel for
            for(omp_int i = 0; i < (omp_int)n; i++)
                for(size_t j = 0; j < m_curr; j++)
                    x0[i] += VT[j][i] * d[j];

            // r
            mul_matrix(x0, r);
            #pragma omp parallel for
            for(omp_int i = 0; i < (omp_int)n; i++)
                r[i] = rp[i] - r[i];
        }
        else
            finished = true;
    }

//    mul_matrix(solution, r);
//    for(size_t i = 0; i < n; i++)
//        r[i] = rp[i] - r[i];
//    discr = sqrt(dot_prod_self(r));
    printf("GMRES_Complex_MKL(%lu) [%d] Residual:\t%5lu\t%.3e\n", (unsigned long)m, numThreads, (unsigned long)iter - 1, discr / rp_norm);

    if(iter >= max_iter)
        printf("Soulution can`t found, iteration limit exceeded!\n");
    else if(residual > eps)
        printf("Soulution can`t found, stagnation detected!\n");
}

GMRES_Complex_MKL::GMRES_Complex_MKL()
{
    d = r = w = x0 = NULL;
    H = VT = NULL;
    m = m_curr = 0;
    ia = ja = NULL;
    aa = NULL;
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
}

GMRES_Complex_MKL::~GMRES_Complex_MKL()
{
    delete [] r;
    delete [] w;
    delete [] d;
    if(H && m)
        for(size_t i = 0; i <= m; i++)
            delete [] H[i];
    delete [] H;
    if(VT && m)
        for(size_t i = 0; i < m; i++)
            delete [] VT[i];
    delete [] VT;
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
