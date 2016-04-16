#include "mkl_stubs.h"
#include <complex>
#include <cmath>
#include <cstdlib>
#include <cstdio>

using namespace std;

namespace mkl_stubs
{

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
    MKL_Complex16 * result = static_cast<MKL_Complex16 *>(dotc);
    result->real = dp_r;
    result->imag = dp_i;
}

void cblas_zdotc_sub(const MKL_INT N, const void * X, const MKL_INT incX, const void * Y, const MKL_INT incY, void * dotc)
{
    const complex<double> * x = static_cast<const complex<double> *>(X);
    const complex<double> * y = static_cast<const complex<double> *>(Y);
    double dp_r = 0.0, dp_i = 0.0;
#pragma omp parallel for reduction(+ : dp_r, dp_i)
    for(MKL_INT i = 0; i < N; i++)
    {
        MKL_INT j = i * incX;
        MKL_INT k = i * incY;
        complex<double> t = conj(x[j]) * y[k];
        dp_r += t.real();
        dp_i += t.imag();
    }
    MKL_Complex16 * result = static_cast<MKL_Complex16 *>(dotc);
    result->real = dp_r;
    result->imag = dp_i;
}

void mkl_zcsrsymv(const char * uplo, const MKL_INT * m, const MKL_Complex16 * a, const MKL_INT * ia,  const MKL_INT * ja, const MKL_Complex16 * x,  MKL_Complex16 * y)
{
    const complex<double> * aa = reinterpret_cast<const complex<double> *>(a);
    const complex<double> * xx = reinterpret_cast<const complex<double> *>(x);
    complex<double> * yy = reinterpret_cast<complex<double> *>(y);
    MKL_INT n = * m;
    if(uplo[0] == 'L' || uplo[0] == 'l')
    {
/*
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
*/
/**/
        complex<double> * mv_tmp;
        MKL_INT * mv_ind, nt1;
        int numThreads = omp_get_max_threads();
        (void)numThreads;

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
/**/
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

int mkl_get_max_threads()
{
    return omp_get_max_threads();
}

}
