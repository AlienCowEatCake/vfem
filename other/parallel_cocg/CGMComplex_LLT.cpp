#include "CGMComplex_LLT.h"
#if defined _MSC_VER && !defined _CRT_SECURE_NO_WARNINGS
#define _CRT_SECURE_NO_WARNINGS
#endif
#include <cstdio>
#include <cmath>

#include <omp.h>
#define THREADS_NUM 4
/*
void CGMComplex_LLT::init(size_t * gi_s, size_t * gj_s, complex<double> * di_s,
                          complex<double> * gg_s, size_t n_s)
{
    gi = gi_s;
    gj = gj_s;
    di = di_s;
    gg = gg_s;
    n = n_s;

    size_t m = gi[n];

    delete [] r;
    delete [] z;
    delete [] p;
    delete [] s;
    delete [] L_di;
    delete [] L_gg;
    delete [] ompvec;

    r = new complex<double> [n];
    z = new complex<double> [n];
    p = new complex<double> [n];
    s = new complex<double> [n];

    L_di = new complex<double> [n];
    L_gg = new complex<double> [m];

    for(size_t i = 0; i < n; i++)
    {
        L_di[i] = di[i];
        //x0[i] = 0.0;
    }

    for(size_t i = 0 ; i < m ; i++)
        L_gg[i] = gg[i];

#if THREADS_NUM > 1
    ompvec = new complex<double> [n * (THREADS_NUM - 1)];
#endif
}

void CGMComplex_LLT::make_LLT_decomposition()
{
    complex<double> sum_d, sum_l;

    for(size_t k = 0; k < n ; k++)
    {
        sum_d = 0;
        size_t i_s = gi[k], i_e = gi[k+1];

        for(size_t i = i_s; i < i_e ; i++)
        {
            sum_l = 0;
            size_t j_s = gi[gj[i]], j_e = gi[gj[i]+1];

            for(size_t m = i_s; m < i; m++)
            {
                for(size_t j = j_s; j < j_e; j++)
                {
                    if(gj[m] == gj[j])
                    {
                        sum_l += L_gg[m] * L_gg[j];
                        j_s++;
                    }
                }
            }
            L_gg[i] = (L_gg[i] -  sum_l) / L_di[gj[i]];

            sum_d += L_gg[i] * L_gg[i];
        }
        L_di[k] = sqrt(L_di[k] - sum_d);
    }
}
*/
complex<double> CGMComplex_LLT::dot_prod_nocj(const complex<double> * a, const complex<double> * b) const
{
    double d_p_re = 0.0, d_p_im = 0.0;
#pragma omp parallel for reduction(+ : d_p_re, d_p_im) num_threads(THREADS_NUM)
    for(int i = 0; i < n; i++) /// WARNING: int в цикле!
    {
        complex<double>d_p(a[i] * b[i]);
        d_p_re += d_p.real();
        d_p_im += d_p.imag();
    }
    return complex<double>(d_p_re, d_p_im);
}

double CGMComplex_LLT::dot_prod_self(const complex<double> * a) const
{
    double d_p = 0.0;
#pragma omp parallel for reduction(+ : d_p) num_threads(THREADS_NUM)
    for(int i = 0; i < n; i++) /// WARNING: int в цикле!
    {
        double re = a[i].real();
        double im = a[i].imag();
        d_p += re * re + im * im;
    }
    return d_p;
}

void CGMComplex_LLT::mul_matrix(const complex<double> * f, complex<double> * x) const
{
#pragma omp parallel num_threads(THREADS_NUM)
{
#pragma omp for
    for(int i = 0; i < n; i++)
        x[i] = di[i] * f[i];

#if THREADS_NUM > 1
#pragma omp for
    for(int i = 0; i < n * (THREADS_NUM - 1); i++)
        ompvec[i] = 0.0;
#endif

#pragma omp for
    for(int i = 0; i < n; i++)
    {
        int rank = omp_get_thread_num();
        if(rank == 0)
        {
            for(size_t k = gi[i], k1 = gi[i+1]; k < k1; k++)
            {
                size_t j = gj[k];
                x[i] += gg[k] * f[j];
                x[j] += gg[k] * f[i];
            }
        }
#if THREADS_NUM > 1
        else
        {
            for(size_t k = gi[i], k1 = gi[i+1]; k < k1; k++)
            {
                size_t j = gj[k];
                size_t adr = (rank - 1) * n;
                ompvec[adr + i] += gg[k] * f[j];
                ompvec[adr + j] += gg[k] * f[i];
            }
        }
#endif
    }

#if THREADS_NUM > 1
#pragma omp for
    for(int i = 0; i < n; i++)
    {
        for(size_t j = 0; j < THREADS_NUM - 1; j++)
        {
            x[i] += ompvec[j * n + i];
        }
    }
#endif
}
}
/*
void CGMComplex_LLT::solve_L(const complex<double> * f, complex<double> * x) const
{
    for(size_t k = 1, k1 = 0; k <= n; k++, k1++)
    {
        complex<double> sum = 0.0;

        for(size_t i = gi[k1]; i < gi[k]; i++)
            sum += L_gg[i] * x[gj[i]];

        x[k1] = (f[k1] - sum) / L_di[k1];
    }
}

void CGMComplex_LLT::solve_LT(complex<double> * f, complex<double> * x) const
{
    for(size_t k = n, k1 = n-1; k > 0; k--, k1--)
    {
        x[k1] = f[k1] / L_di[k1];
        complex<double> v_el = x[k1];

        for(size_t i = gi[k1]; i < gi[k]; i++)
            f[gj[i]] -= L_gg[i] * v_el;
    }
}

void CGMComplex_LLT::solve_LLT(const complex<double> * f, complex<double> * x) const
{
for(size_t i = 0; i < n; i++)
x[i] = f[i];
return;
    solve_L(f, x);
    solve_LT(x, x);
}

bool CGMComplex_LLT::is_fpu_error(double x) const
{
    double y = x - x;
    return x != x || y != y;
}

void CGMComplex_LLT::solve(complex<double> * solution, complex<double> * rp_s, double eps)
{
    // Параметры решателя
    size_t max_iter = 20000;

    rp = rp_s;

    x0 = new complex<double> [n];
    for(size_t i = 0; i < n; i++)
        x0[i] = solution[i];

    mul_matrix(x0, r);
    make_LLT_decomposition();

    for(size_t i = 0; i < n ; i++)
        r[i] = rp[i] - r[i];

    solve_LLT(r, z);
    for(size_t i = 0; i < n; i++)
        p[i] = z[i];

    complex<double> alpha, beta, prod_1, prod_2;
    double discr, rp_norm;

    rp_norm = sqrt(dot_prod_self(rp));
    if(is_fpu_error(rp_norm))
    {
        fprintf(stderr, "Error: FPU error detected in right part!\n");
        for(size_t i = 0; i < n; i++)
            solution[i] = x0[i];
        delete [] x0;
        return;
    }

    prod_1 = dot_prod_nocj(p, r);

    bool finished = false;

    size_t iter;
    for(iter = 1; iter <= max_iter && !finished; iter++)
    {
        discr = sqrt(dot_prod_self(r));
        if(is_fpu_error(discr))
        {
            fprintf(stderr, "Error: FPU error detected in (r, r)!\n");
            for(size_t i = 0; i < n; i++)
                solution[i] = x0[i];
            delete [] x0;
            return;
        }

        double residual = discr / rp_norm;
        if(iter%10 == 0)
        {
            printf("CCGM_LLT Residual:\t%5lu\t%.3e\r", (unsigned long)iter, residual);
            fflush(stdout);
        }

        if(residual > eps)
        {
            mul_matrix(z, s);

            alpha = prod_1 / dot_prod_nocj(s, z);

            for(size_t i = 0; i < n ; i++)
            {
                x0[i] += alpha * z[i];
                r[i] -= alpha * s[i];
            }

            solve_LLT(r, p);
            prod_2 = dot_prod_nocj(p, r);

            beta = prod_2 / prod_1;

            prod_1 = prod_2;

            for(size_t i = 0; i < n; i++)
                z[i] = p[i] + beta * z[i];
        }
        else
            finished = true;
    }

    mul_matrix(x0, r);
    for(size_t i = 0; i < n; i++)
        r[i] = rp[i] - r[i];
    discr = sqrt(dot_prod_self(r));
    printf("CCGM_LLT Residual:\t%5lu\t%.3e\n", (unsigned long)iter, discr / rp_norm);

    if(iter >= max_iter)
        printf("Soulution can`t found, iteration limit exceeded!\n");

    for(size_t i = 0; i < n; i++)
        solution[i] = x0[i];
    delete [] x0;
}

CGMComplex_LLT::CGMComplex_LLT()
{
    r = x0 = z = p = s = L_di = L_gg = ompvec = NULL;
}

CGMComplex_LLT::~CGMComplex_LLT()
{
    delete [] r;
    delete [] z;
    delete [] p;
    delete [] s;
    delete [] L_di;
    delete [] L_gg;
    delete [] ompvec;
}
*/

void CGMComplex_LLT::init(size_t * gi_s, size_t * gj_s, complex<double> * di_s,
                          complex<double> * gg_s, size_t n_s)
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

    r = new complex<double> [n];
    z = new complex<double> [n];
    p = new complex<double> [n];
    s = new complex<double> [n];

ompvec = new complex<double> [n * (THREADS_NUM - 1)];
}
/*
complex<double> CGMComplex_LLT::dot_prod_nocj(const complex<double> * a, const complex<double> * b) const
{
    complex<double> d_p = 0.0;
    for(size_t i = 0; i < n; i++)
        d_p += a[i] * b[i];
    return d_p;
}

double CGMComplex_LLT::dot_prod_self(const complex<double> * a) const
{
    double d_p = 0.0;
    for(size_t i = 0; i < n; i++)
    {
        double re = a[i].real();
        double im = a[i].imag();
        d_p += re * re + im * im;
    }
    return d_p;
}

void CGMComplex_LLT::mul_matrix(const complex<double> * f, complex<double> * x) const
{
    for(size_t i = 0; i < n; i++)
    {
        complex<double> v_el = f[i];
        x[i] = di[i] * v_el;
        for(size_t k = gi[i], k1 = gi[i+1]; k < k1; k++)
        {
            size_t j = gj[k];
            x[i] += gg[k] * f[j];
            x[j] += gg[k] * v_el;
        }
    }
}
*/
bool CGMComplex_LLT::is_fpu_error(double x) const
{
    double y = x - x;
    return x != x || y != y;
}

void CGMComplex_LLT::solve(complex<double> * solution, complex<double> * rp_s, double eps)
{
#if defined SLAE_DEBUG
    FILE * slae_fp = fopen("CGMComplex_LLT_slae_dbg.txt", "w");
#endif

    // Параметры решателя
    size_t max_iter = 15000;

    rp = rp_s;

    mul_matrix(solution, r);
    x0 = new complex<double> [n];
    for(size_t i = 0; i < n; i++)
    {
        x0[i] = solution[i];
        r[i] = rp[i] - r[i];
        z[i] = r[i] / di[i];
        p[i] = z[i];
    }

    complex<double> alpha, beta, prod_1, prod_2;
    double discr, rp_norm;

    rp_norm = sqrt(dot_prod_self(rp));
    if(is_fpu_error(rp_norm))
    {
        fprintf(stderr, "Error: FPU error detected in right part!\n");
        for(size_t i = 0; i < n; i++)
            solution[i] = x0[i];
        delete [] x0;
        return;
    }

    prod_1 = dot_prod_nocj(p, r);

    bool finished = false;

    size_t iter;
    for(iter = 1; iter <= max_iter && !finished; iter++)
    {
        discr = sqrt(dot_prod_self(r));
        if(is_fpu_error(discr))
        {
            fprintf(stderr, "Error: FPU error detected in (r, r)!\n");
            for(size_t i = 0; i < n; i++)
                solution[i] = x0[i];
            delete [] x0;
            return;
        }

        double residual = discr / rp_norm;
        if(iter%10 == 0)
        {
            printf("CGMComplex_LLT Residual:\t%5lu\t%.3e\r", (unsigned long)iter, residual);
            fflush(stdout);
        }
#if defined SLAE_DEBUG
        fprintf(slae_fp, "%lu %.5e\n", (unsigned long)iter, residual);
#endif

        if(residual > eps)
        {
            mul_matrix(z, s);

            alpha = prod_1 / dot_prod_nocj(s, z);

            for(size_t i = 0; i < n ; i++)
            {
                x0[i] += alpha * z[i];
                r[i] -= alpha * s[i];
                p[i] = r[i] / di[i];
            }

            prod_2 = dot_prod_nocj(p, r);

            beta = prod_2 / prod_1;

            prod_1 = prod_2;

            for(size_t i = 0; i < n; i++)
                z[i] = p[i] + beta * z[i];
        }
        else
            finished = true;
    }

    mul_matrix(x0, r);
    for(size_t i = 0; i < n; i++)
        r[i] = rp[i] - r[i];
    discr = sqrt(dot_prod_self(r));
    printf("CGMComplex_LLT Residual:\t%5lu\t%.3e\n", (unsigned long)iter - 1, discr / rp_norm);

    if(iter >= max_iter)
        printf("Soulution can`t found, iteration limit exceeded!\n");

    for(size_t i = 0; i < n; i++)
        solution[i] = x0[i];
    delete [] x0;

#if defined SLAE_DEBUG
    fclose(slae_fp);
#endif
}

CGMComplex_LLT::CGMComplex_LLT()
{
    r = x0 = z = p = s = NULL;
}

CGMComplex_LLT::~CGMComplex_LLT()
{
    delete [] r;
    delete [] z;
    delete [] p;
    delete [] s;
}

