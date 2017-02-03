#if defined(_MSC_VER) && !defined(_CRT_SECURE_NO_WARNINGS)
#define _CRT_SECURE_NO_WARNINGS
#endif
#include "mkl_wrapper.h"
#include "omp_wrapper.h"
#include <complex>
#include <cmath>
#include <cstdlib>
#include <cstdio>

// *************************************************************************************************

#if !defined(USE_MKL)

namespace fem_core { namespace wrappers { namespace mkl_stubs {

void cblas_scopy(const MKL_INT N, const void * X, const MKL_INT incX, void * Y, const MKL_INT incY)
{
    const float * x = static_cast<const float *>(X);
    float * y = static_cast<float *>(Y);
#pragma omp parallel for
    for(MKL_INT i = 0; i < N; i++)
    {
        MKL_INT j = i * incX;
        MKL_INT k = i * incY;
        y[k] = x[j];
    }
}

void cblas_dcopy(const MKL_INT N, const void * X, const MKL_INT incX, void * Y, const MKL_INT incY)
{
    const double * x = static_cast<const double *>(X);
    double * y = static_cast<double *>(Y);
#pragma omp parallel for
    for(MKL_INT i = 0; i < N; i++)
    {
        MKL_INT j = i * incX;
        MKL_INT k = i * incY;
        y[k] = x[j];
    }
}

void cblas_ccopy(const MKL_INT N, const void * X, const MKL_INT incX, void * Y, const MKL_INT incY)
{
    const std::complex<float> * x = static_cast<const std::complex<float> *>(X);
    std::complex<float> * y = static_cast<std::complex<float> *>(Y);
#pragma omp parallel for
    for(MKL_INT i = 0; i < N; i++)
    {
        MKL_INT j = i * incX;
        MKL_INT k = i * incY;
        y[k] = x[j];
    }
}

void cblas_zcopy(const MKL_INT N, const void * X, const MKL_INT incX, void * Y, const MKL_INT incY)
{
    const std::complex<double> * x = static_cast<const std::complex<double> *>(X);
    std::complex<double> * y = static_cast<std::complex<double> *>(Y);
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
    const std::complex<double> * x = static_cast<const std::complex<double> *>(X);
    const std::complex<double> * a = static_cast<const std::complex<double> *>(alpha);
    std::complex<double> * y = static_cast<std::complex<double> *>(Y);
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
    const std::complex<double> * x = static_cast<const std::complex<double> *>(X);
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
    const std::complex<double> * x = static_cast<const std::complex<double> *>(X);
    const std::complex<double> * y = static_cast<const std::complex<double> *>(Y);
    double dp_r = 0.0, dp_i = 0.0;
#pragma omp parallel for reduction(+ : dp_r, dp_i)
    for(MKL_INT i = 0; i < N; i++)
    {
        MKL_INT j = i * incX;
        MKL_INT k = i * incY;
        std::complex<double> t = x[j] * y[k];
        dp_r += t.real();
        dp_i += t.imag();
    }
    MKL_Complex16 * result = static_cast<MKL_Complex16 *>(dotc);
    result->real = dp_r;
    result->imag = dp_i;
}

void cblas_zdotc_sub(const MKL_INT N, const void * X, const MKL_INT incX, const void * Y, const MKL_INT incY, void * dotc)
{
    const std::complex<double> * x = static_cast<const std::complex<double> *>(X);
    const std::complex<double> * y = static_cast<const std::complex<double> *>(Y);
    double dp_r = 0.0, dp_i = 0.0;
#pragma omp parallel for reduction(+ : dp_r, dp_i)
    for(MKL_INT i = 0; i < N; i++)
    {
        MKL_INT j = i * incX;
        MKL_INT k = i * incY;
        std::complex<double> t = conj(x[j]) * y[k];
        dp_r += t.real();
        dp_i += t.imag();
    }
    MKL_Complex16 * result = static_cast<MKL_Complex16 *>(dotc);
    result->real = dp_r;
    result->imag = dp_i;
}

void mkl_zcsrsymv(const char * uplo, const MKL_INT * m, const MKL_Complex16 * a, const MKL_INT * ia,  const MKL_INT * ja, const MKL_Complex16 * x,  MKL_Complex16 * y)
{
    const std::complex<double> * aa = reinterpret_cast<const std::complex<double> *>(a);
    const std::complex<double> * xx = reinterpret_cast<const std::complex<double> *>(x);
    std::complex<double> * yy = reinterpret_cast<std::complex<double> *>(y);
    MKL_INT n = * m;
    if(uplo[0] == 'L' || uplo[0] == 'l')
    {
/*
        for(MKL_INT i = 0; i < n; i++)
        {
            std::complex<double> v_el = xx[i];
            yy[i] = aa[ia[i + 1] - 2] * v_el;
            for(MKL_INT k = ia[i] - 1, k1 = ia[i + 1] - 2; k < k1; k++)
            {
                MKL_INT j = ja[k] - 1;
                yy[i] += aa[k] * xx[j];
                yy[j] += aa[k] * v_el;
            }
        }
*/
/**/
        std::complex<double> * mv_tmp;
        MKL_INT * mv_ind, nt1;
        int numThreads = omp_get_max_threads();
        (void)numThreads;

#pragma omp parallel num_threads(numThreads)
        {
            int numThreads = omp_get_num_threads();
#pragma omp single
            {
                mv_tmp = new std::complex<double> [n * (numThreads - 1)];
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
                yy[i] = std::complex<double>(0.0, 0.0);

#pragma omp for
            for(MKL_INT j = 0; j < (n * nt1); j++)
                mv_tmp[j] = std::complex<double>(0.0, 0.0);

            MKL_INT rank = omp_get_thread_num();
            MKL_INT i_beg = mv_ind[rank];
            MKL_INT i_end = mv_ind[rank + 1];
            for(MKL_INT i = i_beg; i < i_end; i++)
            {
                if(rank == 0)
                {
                    std::complex<double> v_el = xx[i];
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
                    std::complex<double> v_el = xx[i];
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
/**/
        return;
    }
    fprintf(stderr, "mkl_zcsrsymv is not fully implemented!\n");
    exit(-1);
}

void mkl_zcsrtrsv(const char * uplo, const char * transa, const char * diag, const MKL_INT * m, const MKL_Complex16 * a, const MKL_INT * ia, const MKL_INT * ja, const MKL_Complex16 * x, MKL_Complex16 * y)
{
    const std::complex<double> * xx = reinterpret_cast<const std::complex<double> *>(x);
    const std::complex<double> * L_gg = reinterpret_cast<const std::complex<double> *>(a);
    std::complex<double> * yy = reinterpret_cast<std::complex<double> *>(y);
    MKL_INT n = * m;
    if(uplo[0] == 'L' || uplo[0] == 'l' || diag[0] == 'N' || diag[0] == 'n')
    {
        if(transa[0] == 'N' || transa[0] == 'n')
        {
            for(MKL_INT k = 1, k1 = 0; k <= n; k++, k1++)
            {
                std::complex<double> sum = 0.0;
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
                std::complex<double> v_el = yy[k1];
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

int mkl_get_max_threads()
{
    return omp_get_max_threads();
}

}}} // namespace fem_core::wrappers::mkl_stubs

#endif

// *************************************************************************************************

namespace fem_core { namespace wrappers { namespace mkl {

/**
 * @brief Установить количество тредов для MKL равным переменной окружения MKL_NUM_THREADS
 */
void wrapper_mkl_set_env_max_threads()
{
    int num_threads = -1;
    if(num_threads <= 0)
    {
        char * env_num_threads = getenv("MKL_NUM_THREADS");
        if(env_num_threads)
            num_threads = atoi(env_num_threads);
    }
    if(num_threads <= 0)
    {
        char * env_num_threads = getenv("OMP_NUM_THREADS");
        if(env_num_threads)
            num_threads = atoi(env_num_threads);
    }
    if(num_threads <= 0)
    {
        num_threads = 1;
    }
    mkl_set_num_threads(num_threads);
    omp_set_num_threads(num_threads);
}

}}} // namespace fem_core::wrappers::mkl

// *************************************************************************************************

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

// *************************************************************************************************
