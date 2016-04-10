#if defined _MSC_VER && !defined _CRT_SECURE_NO_WARNINGS
#define _CRT_SECURE_NO_WARNINGS
#endif
#include "COCG_GS.h"
#include <cstdio>
#include <cmath>

void COCG_GS::init(const size_t * gi_s, const size_t * gj_s, const complex<double> * di_s,
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

    r = new complex<double> [n];
    z = new complex<double> [n];
    p = new complex<double> [n];
    s = new complex<double> [n];
}

complex<double> COCG_GS::dot_prod_nocj(const complex<double> * a, const complex<double> * b) const
{
    complex<double> d_p = 0.0;
    for(size_t i = 0; i < n; i++)
        d_p += a[i] * b[i];
    return d_p;
}

double COCG_GS::dot_prod_self(const complex<double> * a) const
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

void COCG_GS::mul_matrix(const complex<double> * f, complex<double> * x) const
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

void COCG_GS::solve_GS(const complex<double> * f, complex<double> * x) const
{
    double w = 1.0;

    complex<double> * tmp = new complex<double> [n];
    for(size_t i = 0; i < n; i++)
        tmp[i] = f[i];

    // (D + L)^-1 tmp = x
    for(size_t k = 1, k1 = 0; k <= n; k++, k1++)
    {
        complex<double> sum = 0.0;
        for(size_t i = gi[k1]; i < gi[k]; i++)
            sum += gg[i] * x[gj[i]] * w;
        x[k1] = (tmp[k1] - sum) / di[k1];
    }

    // tmp = D * x
    for(size_t i = 0; i < n; i++)
        tmp[i] = x[i] * di[i];

    // (D + LT)^-1 * tmp = x
    for(size_t k = n, k1 = n-1; k > 0; k--, k1--)
    {
        x[k1] = tmp[k1] / di[k1];
        complex<double> v_el = x[k1];
        for(size_t i = gi[k1]; i < gi[k]; i++)
            tmp[gj[i]] -= gg[i] * v_el * w;
    }

    delete [] tmp;
}

bool COCG_GS::is_fpu_error(double x) const
{
    double y = x - x;
    return x != x || y != y;
}

void COCG_GS::solve(complex<double> * solution, const complex<double> * rp_s, double eps, size_t max_iter)
{
    rp = rp_s;

    x0 = new complex<double> [n];
    for(size_t i = 0; i < n; i++)
        x0[i] = solution[i];

    mul_matrix(x0, r);

    for(size_t i = 0; i < n ; i++)
        r[i] = rp[i] - r[i];

    solve_GS(r, z);
    for(size_t i = 0; i < n; i++)
        p[i] = z[i];

    complex<double> alpha, beta, prod_1, prod_2;
    double discr = 0.0, rp_norm;

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
    for(iter = 0; iter <= max_iter && !finished; iter++)
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
            printf("COCG_GS Residual:\t%5lu\t%.3e\r", (unsigned long)iter, residual);
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

            solve_GS(r, p);
            prod_2 = dot_prod_nocj(p, r);

            beta = prod_2 / prod_1;

            prod_1 = prod_2;

            for(size_t i = 0; i < n; i++)
                z[i] = p[i] + beta * z[i];
        }
        else
            finished = true;
    }

//    mul_matrix(x0, r);
//    for(size_t i = 0; i < n; i++)
//        r[i] = rp[i] - r[i];
//    discr = sqrt(dot_prod_self(r));
    printf("COCG_GS Residual:\t%5lu\t%.3e\n", (unsigned long)iter - 1, discr / rp_norm);

    if(iter >= max_iter)
        printf("Soulution can`t found, iteration limit exceeded!\n");

    for(size_t i = 0; i < n; i++)
        solution[i] = x0[i];
    delete [] x0;
}

COCG_GS::COCG_GS()
{
    r = x0 = z = p = s = NULL;
}

COCG_GS::~COCG_GS()
{
    delete [] r;
    delete [] z;
    delete [] p;
    delete [] s;
}
