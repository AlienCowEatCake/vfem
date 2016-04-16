#if defined(_MSC_VER) && !defined(_CRT_SECURE_NO_WARNINGS)
#define _CRT_SECURE_NO_WARNINGS
#endif
#include "COCG_Di_Smooth_OpenMP.h"
#include <cstdio>
#include <cmath>

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

void COCG_Di_Smooth_OpenMP::init(const size_t * gi_s, const size_t * gj_s, const complex<double> * di_s,
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
    delete [] mv_tmp;
    delete [] mv_ind;

    r = new complex<double> [n];
    z = new complex<double> [n];
    p = new complex<double> [n];
    s = new complex<double> [n];
    xs = new complex<double> [n];
    rs = new complex<double> [n];
    mv_tmp = new complex<double> [n * (size_t)(numThreads - 1)];

    // Индексы начала/конца в f и x
    mv_ind = new size_t [numThreads + 1];
    // Начинает первый в начале
    mv_ind[0] = 0;
    // Заканчивает последний в конце
    mv_ind[numThreads] = n;
    // Выравнивание по строкам
    for(size_t i = 0, curr_ind = 1, curr_size = gi[n] / (size_t)numThreads; i < n; i++)
    {
        // Если где-то между текущей и следующей строкой проходит разбиение
        if(gi[i + 1] >= curr_size && gi[i] <= curr_size)
        {
            // Найдем куда ближе сдвигать - к текущей или к следующей
            // Если ближе к следующей - сделаем ее текущей
            if(gi[i + 1] - curr_size < curr_size - gi[i])
                i++;
            // Сдвигаемся
            mv_ind[curr_ind] = i;
            curr_size += gi[n] / (size_t)numThreads;
            curr_ind++;
        }
    }
}

complex<double> COCG_Di_Smooth_OpenMP::dot_prod_nocj(const complex<double> * a, const complex<double> * b) const
{
    double d_p_real = 0.0;
    double d_p_imag = 0.0;
#pragma omp parallel for reduction(+ : d_p_real, d_p_imag)
    for(omp_int i = 0; i < (omp_int)n; i++)
    {
        complex<double> c = a[i] * b[i];
        d_p_real += c.real();
        d_p_imag += c.imag();
    }
    return complex<double>(d_p_real, d_p_imag);
}

double COCG_Di_Smooth_OpenMP::dot_prod_self(const complex<double> * a) const
{
    double d_p = 0.0;
#pragma omp parallel for reduction(+ : d_p)
    for(omp_int i = 0; i < (omp_int)n; i++)
    {
        double re = a[i].real();
        double im = a[i].imag();
        d_p += re * re + im * im;
    }
    return d_p;
}

double COCG_Di_Smooth_OpenMP::dot_prod_real(const complex<double> * a, const complex<double> * b) const
{
    double d_p = 0.0;
#pragma omp parallel for reduction(+ : d_p)
    for(omp_int i = 0; i < (omp_int)n; i++)
        d_p += a[i].real() * b[i].real() + a[i].imag() * b[i].imag();
    return d_p;
}

void COCG_Di_Smooth_OpenMP::mul_matrix(const complex<double> * f, complex<double> * x) const
{
    size_t nt1 = (size_t)(numThreads - 1);

#pragma omp parallel num_threads(numThreads)
    {
#pragma omp for nowait
        for(omp_int i = 0; i < (omp_int)n; i++)
            x[i] = complex<double>(0.0, 0.0);

#pragma omp for
        for(omp_int j = 0; j < (omp_int)(n * nt1); j++)
            mv_tmp[j] = complex<double>(0.0, 0.0);

        omp_int rank = (omp_int)omp_get_thread_num();
        omp_int i_beg = (omp_int)mv_ind[rank];
        omp_int i_end = (omp_int)mv_ind[rank + 1];
        for(omp_int i = i_beg; i < i_end; i++)
        {
            if(rank == 0)
            {
                complex<double> v_el = f[i];
                x[i] = di[i] * v_el;
                for(size_t k = gi[i], k1 = gi[i + 1]; k < k1; k++)
                {
                    size_t j = gj[k];
                    x[i] += gg[k] * f[j];
                    x[j] += gg[k] * v_el;
                }
            }
            else
            {
                size_t adr = (rank - 1) * n;
                size_t adr_i = adr + i;
                complex<double> v_el = f[i];
                mv_tmp[adr_i] = di[i] * v_el;
                for(size_t k = gi[i], k1 = gi[i + 1]; k < k1; k++)
                {
                    size_t j = gj[k];
                    mv_tmp[adr_i] += gg[k] * f[j];
                    mv_tmp[adr + j] += gg[k] * v_el;
                }
            }
        }
#pragma omp barrier

#pragma omp for
        for(omp_int i = 0; i < (omp_int)n; i++)
        {
            for(size_t j = 0; j < nt1; j++)
                x[i] += mv_tmp[j * n + i];
        }
    }
}

bool COCG_Di_Smooth_OpenMP::is_fpu_error(double x) const
{
    double y = x - x;
    return x != x || y != y;
}

void COCG_Di_Smooth_OpenMP::solve(complex<double> * solution, const complex<double> * rp_s,
                                  double eps, size_t max_iter)
{
    eps *= eps;

    rp = rp_s;

    mul_matrix(solution, r);
    x0 = new complex<double> [n];

#pragma omp parallel for
    for(omp_int i = 0; i < (omp_int)n; i++)
    {
        xs[i] = x0[i] = solution[i];
        rs[i] = r[i] = rp[i] - r[i];
        z[i] = r[i] / di[i];
        p[i] = z[i];
    }

    complex<double> alpha, beta, prod_1, prod_2;
    double discr = 0.0, rp_norm, eta;

    rp_norm = dot_prod_self(rp);
    if(is_fpu_error(rp_norm))
    {
        fprintf(stderr, "Error: FPU error detected in right part!\n");
#pragma omp parallel for
        for(omp_int i = 0; i < (omp_int)n; i++)
            solution[i] = x0[i];
        delete [] x0;
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
#pragma omp parallel for
            for(omp_int i = 0; i < (omp_int)n; i++)
                solution[i] = xs[i];
            delete [] x0;
            return;
        }

        double residual = discr / rp_norm;
        if(iter%10 == 0)
        {
            printf("COCG_Di_Smooth_OpenMP [%d] Residual:\t%5lu\t%.3e\r", numThreads, (unsigned long)iter, sqrt(residual));
            fflush(stdout);
        }

        if(residual > eps)
        {
            mul_matrix(z, s);

            alpha = prod_1 / dot_prod_nocj(s, z);

#pragma omp parallel for
            for(omp_int i = 0; i < (omp_int)n ; i++)
            {
                x0[i] += alpha * z[i];
                r[i] -= alpha * s[i];
                p[i] = r[i] / di[i];
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

#pragma omp parallel for
            for(omp_int i = 0; i < (omp_int)n; i++)
                z[i] = p[i] + beta * z[i];
        }
        else
            finished = true;
    }

//    mul_matrix(xs, r);
//#pragma omp parallel for
//    for(omp_int i = 0; i < (omp_int)n; i++)
//        r[i] = rp[i] - r[i];
//    discr = dot_prod_self(r);
    printf("COCG_Di_Smooth_OpenMP [%d] Residual:\t%5lu\t%.3e\n", numThreads, (unsigned long)iter - 1, sqrt(discr / rp_norm));

    if(iter >= max_iter)
        printf("Soulution can`t found, iteration limit exceeded!\n");

#pragma omp parallel for
    for(omp_int i = 0; i < (omp_int)n; i++)
        solution[i] = xs[i];
    delete [] x0;
}

COCG_Di_Smooth_OpenMP::COCG_Di_Smooth_OpenMP()
{
    numThreads = -1;
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
    omp_set_num_threads(numThreads);
    numThreads = omp_get_max_threads();
    r = x0 = z = p = s = xs = rs = NULL;
    mv_tmp = NULL;
    mv_ind = NULL;
}

COCG_Di_Smooth_OpenMP::~COCG_Di_Smooth_OpenMP()
{
    delete [] r;
    delete [] z;
    delete [] p;
    delete [] s;
    delete [] xs;
    delete [] rs;
    delete [] mv_tmp;
    delete [] mv_ind;
}
