#if defined(_MSC_VER) && !defined(_CRT_SECURE_NO_WARNINGS)
#define _CRT_SECURE_NO_WARNINGS
#endif
#include "COCR_Di.h"
#include <cstdio>
#include <cmath>

void COCR_Di::init(const size_t * gi_s, const size_t * gj_s, const complex<double> * di_s,
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
    delete [] w;
    delete [] a;

    r = new complex<double> [n];
    z = new complex<double> [n];
    p = new complex<double> [n];
    s = new complex<double> [n];
    w = new complex<double> [n];
    a = new complex<double> [n];
}

complex<double> COCR_Di::dot_prod_nocj(const complex<double> * a, const complex<double> * b) const
{
    complex<double> d_p = 0.0;
    for(size_t i = 0; i < n; i++)
        d_p += a[i] * b[i];
    return d_p;
}

double COCR_Di::dot_prod_self(const complex<double> * a) const
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

void COCR_Di::mul_matrix(const complex<double> * f, complex<double> * x) const
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

bool COCR_Di::is_fpu_error(double x) const
{
    double y = x - x;
    return x != x || y != y;
}

void COCR_Di::solve(complex<double> * solution, const complex<double> * rp_s,
                    double eps, size_t max_iter)
{
    rp = rp_s;

    x0 = new complex<double> [n];
    for(size_t i = 0; i < n; i++)
        x0[i] = solution[i];

    mul_matrix(x0, r);

    for(size_t i = 0; i < n ; i++)
        r[i] = rp[i] - r[i];

    for(size_t i = 0; i < n; i++) // solve_D
        s[i] = r[i] / di[i];
    mul_matrix(s, z);
    for(size_t i = 0; i < n; i++) // solve_D
        w[i] = z[i] / di[i];
    for(size_t i = 0; i < n; i++)
    {
        p[i] = s[i];
        a[i] = z[i];
    }
    complex<double> dp_as = dot_prod_nocj(a, s), dp_as_new;
    complex<double> alpha, beta;
    double discr = 0.0, rp_norm, residual = 0.0, residual_old;

    rp_norm = sqrt(dot_prod_self(rp));
    if(is_fpu_error(rp_norm))
    {
        fprintf(stderr, "Error: FPU error detected in right part!\n");
        for(size_t i = 0; i < n; i++)
            solution[i] = x0[i];
        delete [] x0;
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
            for(size_t i = 0; i < n; i++)
                solution[i] = x0[i];
            delete [] x0;
            return;
        }

        residual_old = residual;
        residual = discr / rp_norm;
        if(iter%10 == 0)
        {
            printf("COCR_Di Residual:\t%5lu\t%.3e\r", (unsigned long)iter, residual);
            fflush(stdout);
        }

        if(residual > eps && fabs(residual - residual_old) / (residual) > 1e-15)
        {
            alpha = dp_as / dot_prod_nocj(z, w);
            for(size_t i = 0; i < n ; i++)
            {
                x0[i] += alpha * p[i];
                r[i] -= alpha * z[i];
                s[i] -= alpha * w[i];
            }
            mul_matrix(s, a);
            dp_as_new = dot_prod_nocj(a, s);
            beta = dp_as_new / dp_as;
            dp_as = dp_as_new;
            for(size_t i = 0; i < n ; i++)
            {
                p[i] = s[i] + beta * p[i];
                z[i] = a[i] + beta * z[i];
                w[i] = z[i] / di[i];    // solve_D
            }
        }
        else
            finished = true;
    }

//    mul_matrix(x0, r);
//    for(size_t i = 0; i < n; i++)
//        r[i] = rp[i] - r[i];
//    discr = sqrt(dot_prod_self(r));
    printf("COCR_Di Residual:\t%5lu\t%.3e\n", (unsigned long)iter - 1, discr / rp_norm);

    if(iter >= max_iter)
        printf("Soulution can`t found, iteration limit exceeded!\n");
    else if(residual > eps)
        printf("Soulution can`t found, stagnation detected!\n");

    for(size_t i = 0; i < n; i++)
        solution[i] = x0[i];
    delete [] x0;
}

COCR_Di::COCR_Di()
{
    r = x0 = z = p = s = w = a = NULL;
}

COCR_Di::~COCR_Di()
{
    delete [] r;
    delete [] z;
    delete [] p;
    delete [] s;
    delete [] w;
    delete [] a;
}
