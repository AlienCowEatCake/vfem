#if defined _MSC_VER && !defined _CRT_SECURE_NO_WARNINGS
#define _CRT_SECURE_NO_WARNINGS
#endif
#include "COCG_LLT_Smooth_MKL.h"
#include <cstdio>
#include <cmath>
#include <algorithm>

#define PRECONDITIONER_NONE 0x01
#define PRECONDITIONER_DI   0x02
#define PRECONDITIONER_LLT  0x03

#if !defined PRECONDITIONER
#define PRECONDITIONER PRECONDITIONER_LLT
#endif

#if defined _MSC_VER
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

void COCG_LLT_Smooth_MKL::init(const size_t * gi_s, const size_t * gj_s, const complex<double> * di_s,
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

#if PRECONDITIONER == PRECONDITIONER_LLT
    delete [] L_aa;
    delete [] LLT_tmp;
    LLT_tmp = new MKL_Complex16 [n];
    L_aa = new MKL_Complex16 [gi[n] + n];
    make_LLT_decomposition();
#endif
}

void COCG_LLT_Smooth_MKL::make_LLT_decomposition()
{
    cblas_zcopy(ia[m] - 1, aa, 1, L_aa, 1);
    complex<double> * L_gg = reinterpret_cast<complex<double> *>(L_aa);

    complex<double> sum_d, sum_l;

    for(size_t k = 0; k < n ; k++)
    {
        sum_d = 0;
        size_t i_s = ia[k] - 1, i_e = ia[k+1] - 2;
        for(size_t i = i_s; i < i_e ; i++)
        {
            sum_l = 0;
            size_t j_s = ia[ja[i]-1]-1, j_e = ia[ja[i]]-2;

            for(size_t m = i_s; m < i; m++)
            {
                for(size_t j = j_s; j < j_e; j++)
                {
                    if(ja[m] == ja[j])
                    {
                        sum_l += L_gg[m] * L_gg[j];
                        j_s++;
                    }
                }
            }
            L_gg[i] = (L_gg[i] -  sum_l) / L_gg[j_e];

            sum_d += L_gg[i] * L_gg[i];
        }
        size_t i_di = ia[k+1]-2;
        L_gg[i_di] = sqrt(L_gg[i_di] - sum_d);
    }
}

complex<double> COCG_LLT_Smooth_MKL::dot_prod_nocj(const complex<double> * a, const complex<double> * b) const
{
    complex<double> result;
    MKL_Complex16 * a1 = reinterpret_cast<MKL_Complex16 *>(const_cast<complex<double> *>(a));
    MKL_Complex16 * b1 = reinterpret_cast<MKL_Complex16 *>(const_cast<complex<double> *>(b));
    MKL_Complex16 * pres = reinterpret_cast<MKL_Complex16 *>(&result);
    cblas_zdotu_sub(m, a1, 1, b1, 1, pres);
    return result;
}

double COCG_LLT_Smooth_MKL::dot_prod_self(const complex<double> * a) const
{
    MKL_Complex16 * a1 = reinterpret_cast<MKL_Complex16 *>(const_cast<complex<double> *>(a));
    double result = cblas_dznrm2(m, a1, 1);
    return result * result;
}

double COCG_LLT_Smooth_MKL::dot_prod_real(const complex<double> * a, const complex<double> * b) const
{
    double d_p = 0.0;
#pragma omp parallel for reduction(+ : d_p)
    for(omp_int i = 0; i < (omp_int)n; i++)
        d_p += a[i].real() * b[i].real() + a[i].imag() * b[i].imag();
    return d_p;
}

void COCG_LLT_Smooth_MKL::mul_matrix(const complex<double> * f, complex<double> * x) const
{
    MKL_Complex16 * in_v = reinterpret_cast<MKL_Complex16 *>(const_cast<complex<double> *>(f));
    MKL_Complex16 * out_v = reinterpret_cast<MKL_Complex16 *>(const_cast<complex<double> *>(x));
    mkl_zcsrsymv("L", &m, aa, ia, ja, in_v, out_v);
}

void COCG_LLT_Smooth_MKL::solve_LLT(const complex<double> * f, complex<double> * x) const
{
    const MKL_Complex16 * ff = reinterpret_cast<const MKL_Complex16 *>(f);
    MKL_Complex16 * xx = reinterpret_cast<MKL_Complex16 *>(x);
    mkl_zcsrtrsv("L", "N", "N", &m, L_aa, ia, ja, ff, LLT_tmp);
    mkl_zcsrtrsv("L", "T", "N", &m, L_aa, ia, ja, LLT_tmp, xx);
}

bool COCG_LLT_Smooth_MKL::is_fpu_error(double x) const
{
    double y = x - x;
    return x != x || y != y;
}

void COCG_LLT_Smooth_MKL::solve(complex<double> * solution, const complex<double> * rp_s, double eps, size_t max_iter)
{
    eps *= eps;

    rp = rp_s;

    x0 = solution;
    mul_matrix(x0, r);

#pragma omp parallel for
    for(omp_int i = 0; i < (omp_int)n; i++)
    {
        r[i] = rp[i] - r[i];
#if PRECONDITIONER == PRECONDITIONER_DI
        z[i] = r[i] / di[i];
#endif
#if PRECONDITIONER == PRECONDITIONER_NONE
        z[i] = r[i];
#endif
    }
#if PRECONDITIONER == PRECONDITIONER_LLT
    solve_LLT(r, z);
#endif
    normal_zcopy(n, solution, xs);
    normal_zcopy(n, r, rs);
    normal_zcopy(n, z, p);

    complex<double> alpha, beta, prod_1, prod_2;
    double discr, rp_norm, eta;

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
#if PRECONDITIONER == PRECONDITIONER_LLT
            printf("COCG_LLT_Smooth_MKL [%d] Residual:\t%5lu\t%.3e\r", numThreads, (unsigned long)iter, sqrt(residual));
#endif
#if PRECONDITIONER == PRECONDITIONER_DI
            printf("COCG_Di_Smooth_MKL [%d] Residual:\t%5lu\t%.3e\r", numThreads, (unsigned long)iter, sqrt(residual));
#endif
#if PRECONDITIONER == PRECONDITIONER_NONE
            printf("COCG_Smooth_MKL [%d] Residual:\t%5lu\t%.3e\r", numThreads, (unsigned long)iter, sqrt(residual));
#endif

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
#if PRECONDITIONER == PRECONDITIONER_DI
                p[i] = r[i] / di[i];
#endif
#if PRECONDITIONER == PRECONDITIONER_NONE
                p[i] = r[i];
#endif
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

#if PRECONDITIONER == PRECONDITIONER_LLT
            solve_LLT(r, p);
#endif
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
#if PRECONDITIONER == PRECONDITIONER_LLT
    printf("COCG_LLT_Smooth_MKL [%d] Residual:\t%5lu\t%.3e\n", numThreads, (unsigned long)iter - 1, sqrt(discr / rp_norm));
#endif
#if PRECONDITIONER == PRECONDITIONER_DI
    printf("COCG_Di_Smooth_MKL [%d] Residual:\t%5lu\t%.3e\n", numThreads, (unsigned long)iter - 1, sqrt(discr / rp_norm));
#endif
#if PRECONDITIONER == PRECONDITIONER_NONE
    printf("COCG_Smooth_MKL [%d] Residual:\t%5lu\t%.3e\n", numThreads, (unsigned long)iter - 1, sqrt(discr / rp_norm));
#endif

    if(iter >= max_iter)
        printf("Soulution can`t found, iteration limit exceeded!\n");

    normal_zcopy(n, xs, solution);
}

COCG_LLT_Smooth_MKL::COCG_LLT_Smooth_MKL()
{
    numThreads = -1;
    if(numThreads <= 0)
    {
        char * env_num_threads;
        env_num_threads = getenv("MKL_NUM_THREADS");
        if(env_num_threads)
            sscanf(env_num_threads, "%d", &numThreads);
    }
    if(numThreads <= 0)
    {
        char * env_num_threads;
        env_num_threads = getenv("OMP_NUM_THREADS");
        if(env_num_threads)
            sscanf(env_num_threads, "%d", &numThreads);
    }
    if(numThreads <= 0)
    {
        numThreads = 1;
    }
    mkl_set_num_threads(numThreads);
    omp_set_num_threads(numThreads);
    r = x0 = z = p = s = xs = rs = NULL;
    L_aa = LLT_tmp = NULL;
    ia = ja = NULL;
    aa = NULL;
}

COCG_LLT_Smooth_MKL::~COCG_LLT_Smooth_MKL()
{
    delete [] r;
    delete [] z;
    delete [] p;
    delete [] s;
    delete [] L_aa;
    delete [] LLT_tmp;
    delete [] xs;
    delete [] rs;
    delete [] ia;
    delete [] ja;
    delete [] aa;
}

// =================================================================================================

#if !defined USE_MKL

void cblas_zcopy(const MKL_INT N, const void * X, const MKL_INT incX, void * Y, const MKL_INT incY)
{
    const complex<double> * x = static_cast<const complex<double> *>(X);
    complex<double> * y = static_cast<complex<double> *>(Y);
#pragma omp parallel for
    for(MKL_INT i = 0; i < N; i++)
    {
        MKL_INT j = i * incX;
        MKL_INT k = i * incY;
        y[k] = x[j];
    }
}

void cblas_zaxpy(const MKL_INT N, const void * alpha, const void * X, const MKL_INT incX, void * Y, const MKL_INT incY)
{
    const complex<double> * x = static_cast<const complex<double> *>(X);
    const complex<double> * a = static_cast<const complex<double> *>(alpha);
    complex<double> * y = static_cast<complex<double> *>(Y);
#pragma omp parallel for
    for(MKL_INT i = 0; i < N; i++)
    {
        MKL_INT j = i * incX;
        MKL_INT k = i * incY;
        y[k] += * a * x[j];
    }
}

double cblas_dznrm2(const MKL_INT N, const void * X, const MKL_INT incX)
{
    const complex<double> * x = static_cast<const complex<double> *>(X);
    double result = 0.0;
#pragma omp parallel for reduction(+ : result)
    for(MKL_INT i = 0; i < N; i++)
    {
        MKL_INT j = i * incX;
        double re = x[j].real();
        double im = x[j].imag();
        result += re * re + im * im;
    }
    return sqrt(result);
}

void cblas_zdotu_sub(const MKL_INT N, const void * X, const MKL_INT incX, const void * Y, const MKL_INT incY, void * dotc)
{
    const complex<double> * x = static_cast<const complex<double> *>(X);
    const complex<double> * y = static_cast<const complex<double> *>(Y);
    double dp_r = 0.0, dp_i = 0.0;
#pragma omp parallel for reduction(+ : dp_r, dp_i)
    for(MKL_INT i = 0; i < N; i++)
    {
        MKL_INT j = i * incX;
        MKL_INT k = i * incY;
        complex<double> t = x[j] * y[k];
        dp_r += t.real();
        dp_i += t.imag();
    }
    complex<double> * result = static_cast<complex<double> *>(dotc);
    result->real(dp_r);
    result->imag(dp_i);
}

void mkl_zcsrsymv(const char * uplo, const MKL_INT * m, const MKL_Complex16 * a, const MKL_INT * ia,  const MKL_INT * ja, const MKL_Complex16 * x,  MKL_Complex16 * y)
{
    const complex<double> * aa = reinterpret_cast<const complex<double> *>(a);
    const complex<double> * xx = reinterpret_cast<const complex<double> *>(x);
    complex<double> * yy = reinterpret_cast<complex<double> *>(y);
    MKL_INT n = * m;
    if(uplo[0] == 'L' || uplo[0] == 'l')
    {
/**/
        for(MKL_INT i = 0; i < n; i++)
        {
            complex<double> v_el = xx[i];
            yy[i] = aa[ia[i + 1] - 2] * v_el;
            for(MKL_INT k = ia[i] - 1, k1 = ia[i + 1] - 2; k < k1; k++)
            {
                MKL_INT j = ja[k] - 1;
                yy[i] += aa[k] * xx[j];
                yy[j] += aa[k] * v_el;
            }
        }
/**/
/*
        complex<double> * mv_tmp;
        MKL_INT * mv_ind, nt1;
        int numThreads = omp_get_num_threads();

#pragma omp parallel num_threads(numThreads)
        {
            int numThreads = omp_get_num_threads();
#pragma omp single
            {
                mv_tmp = new complex<double> [n * (numThreads - 1)];
                mv_ind = new MKL_INT [numThreads + 1];
                mv_ind[0] = 0;
                mv_ind[numThreads] = n;
                for(MKL_INT i = 0, curr_ind = 1, curr_size = ia[n] / numThreads; i < n; i++)
                {
                    if(ia[i + 1] >= curr_size && ia[i] <= curr_size)
                    {
                        if(ia[i + 1] - curr_size < curr_size - ia[i]) i++;
                        mv_ind[curr_ind] = i;
                        curr_size += ia[n] / numThreads;
                        curr_ind++;
                    }
                }
                nt1 = numThreads - 1;
            }

#pragma omp for nowait
            for(MKL_INT i = 0; i < n; i++)
                yy[i] = complex<double>(0.0, 0.0);

#pragma omp for
            for(MKL_INT j = 0; j < (n * nt1); j++)
                mv_tmp[j] = complex<double>(0.0, 0.0);

            MKL_INT rank = omp_get_thread_num();
            MKL_INT i_beg = mv_ind[rank];
            MKL_INT i_end = mv_ind[rank + 1];
            for(MKL_INT i = i_beg; i < i_end; i++)
            {
                if(rank == 0)
                {
                    complex<double> v_el = xx[i];
                    yy[i] = aa[ia[i + 1] - 2] * v_el;
                    for(MKL_INT k = ia[i] - 1, k1 = ia[i + 1] - 2; k < k1; k++)
                    {
                        MKL_INT j = ja[k] - 1;
                        yy[i] += aa[k] * xx[j];
                        yy[j] += aa[k] * v_el;
                    }
                }
                else
                {
                    MKL_INT adr = (rank - 1) * n;
                    MKL_INT adr_i = adr + i;
                    complex<double> v_el = xx[i];
                    mv_tmp[adr_i] = aa[ia[i + 1] - 2] * v_el;
                    for(MKL_INT k = ia[i] - 1, k1 = ia[i + 1] - 2; k < k1; k++)
                    {
                        MKL_INT j = ja[k] - 1;
                        mv_tmp[adr_i] += aa[k] * xx[j];
                        mv_tmp[adr + j] += aa[k] * v_el;
                    }
                }
            }
#pragma omp barrier

#pragma omp for
            for(MKL_INT i = 0; i < n; i++)
            {
                for(MKL_INT j = 0; j < nt1; j++)
                    yy[i] += mv_tmp[j * n + i];
            }
        }

        delete [] mv_tmp;
        delete [] mv_ind;
*/
        return;
    }
    fprintf(stderr, "mkl_zcsrsymv is not fully implemented!\n");
    exit(-1);
}

void mkl_zcsrtrsv(const char * uplo, const char * transa, const char * diag, const MKL_INT * m, const MKL_Complex16 * a, const MKL_INT * ia, const MKL_INT * ja, const MKL_Complex16 * x, MKL_Complex16 * y)
{
    const complex<double> * xx = reinterpret_cast<const complex<double> *>(x);
    const complex<double> * L_gg = reinterpret_cast<const complex<double> *>(a);
    complex<double> * yy = reinterpret_cast<complex<double> *>(y);
    MKL_INT n = * m;
    if(uplo[0] == 'L' || uplo[0] == 'l' || diag[0] == 'N' || diag[0] == 'n')
    {
        if(transa[0] == 'N' || transa[0] == 'n')
        {
            for(MKL_INT k = 1, k1 = 0; k <= n; k++, k1++)
            {
                complex<double> sum = 0.0;
                MKL_INT i_b = ia[k1]-1, i_e = ia[k]-2;
                for(MKL_INT i = i_b; i < i_e; i++)
                    sum += L_gg[i] * yy[ja[i]-1];
                MKL_INT i_d = ia[k]-2;
                yy[k1] = (xx[k1] - sum) / L_gg[i_d];
            }
            return;
        }
        else if(transa[0] == 'T' || transa[0] == 't' || transa[0] == 'C' || transa[0] == 'c')
        {
            for(MKL_INT k = 0; k < n; k++)
                yy[k] = xx[k];
            for(MKL_INT k = n, k1 = n-1; k > 0; k--, k1--)
            {
                MKL_INT i_d = ia[k]-2;
                yy[k1] = yy[k1] / L_gg[i_d];
                complex<double> v_el = yy[k1];
                MKL_INT i_b = ia[k1]-1, i_e = ia[k]-2;
                for(MKL_INT i = i_b; i < i_e; i++)
                    yy[ja[i]-1] -= L_gg[i] * v_el;
            }
            return;
        }
    }
    fprintf(stderr, "mkl_zcsrtrsv is not fully implemented!\n");
    exit(-1);
}

void mkl_set_num_threads(int n)
{
    omp_set_num_threads(n);
}

#endif

// =================================================================================================

#if defined USE_MKL && defined _MSC_VER

#define MKL_TBB_THREADS
//#define MKL_OMP_THREADS

#if defined _WIN32
#if !defined _WIN64
#pragma comment(lib, "mkl_intel_c.lib")
#else
#pragma comment(lib, "mkl_intel_lp64.lib")
#endif
#pragma comment(lib, "mkl_core.lib")
#if defined MKL_TBB_THREADS
#pragma comment(lib, "mkl_tbb_thread.lib")
#pragma comment(lib, "tbb.lib")
#elif defined MKL_OMP_THREADS
#pragma comment(lib, "mkl_intel_thread.lib")
#pragma comment(lib, "libiomp5md.lib")
#else
#pragma comment(lib, "mkl_sequential.lib")
#endif
#endif

#endif
