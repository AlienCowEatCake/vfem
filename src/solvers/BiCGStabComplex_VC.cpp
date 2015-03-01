#include "BiCGStabComplex_VC.h"
#if defined _MSC_VER && !defined _CRT_SECURE_NO_WARNINGS
#define _CRT_SECURE_NO_WARNINGS
#endif
#include <cstdio>
#include <cmath>

BiCGStabComplex_VC::BiCGStabComplex_VC()
{
    r = v = s = p = t = r2 = NULL;
}

BiCGStabComplex_VC::~BiCGStabComplex_VC()
{
    delete [] r;
    delete [] v;
    delete [] s;
    delete [] p;
    delete [] t;
    delete [] r2;
}

void BiCGStabComplex_VC::init(size_t * gi_s, size_t * gj_s, complex<double> * di_s,
                              complex<double> * gg_s, size_t n_s)
{
    gi = gi_s;
    gj = gj_s;
    di = di_s;
    gg = gg_s;
    n = n_s;

    delete [] r;
    delete [] v;
    delete [] s;
    delete [] p;
    delete [] t;
    delete [] r2;

    r = new complex<double> [n];
    v = new complex<double> [n];
    s = new complex<double> [n];
    p = new complex<double> [n];
    t = new complex<double> [n];
    r2 = new complex<double> [n];
}

complex<double> BiCGStabComplex_VC::dot_prod(const complex<double> * a, const complex<double> * b) const
{
    complex<double> res;
    for(size_t i = 0; i < n; i++)
    {
        res += conj(a[i]) * b[i];
    }
    return res;
}

complex<double> BiCGStabComplex_VC::dot_prod_nocj(const complex<double> * a, const complex<double> * b) const
{
    complex<double> res;
    for(size_t i = 0; i < n; i++)
    {
        res += a[i] * b[i];
    }
    return res;
}

void BiCGStabComplex_VC::mul_matrix(const complex<double> * f, complex<double> * x) const
{
    for(size_t i = 0; i < n; i++)
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
}

void BiCGStabComplex_VC::solve(complex<double> *solution, complex<double> *rp, double gamma)
{
    double eps = gamma;
    size_t max_iter = /*(size_t)sqrt(n)*/ 15000;

    complex<double> omega, alpha, ro, ro_prev, beta;
    double rp_norm = sqrt(dot_prod(rp, rp).real());
    x0 = new complex<double> [n];

    mul_matrix(solution, t);
    for(size_t i = 0; i < n; i++)
    {
        x0[i] = solution[i];
        r[i] = rp[i] - t[i];
        r2[i] = 1.0;
        v[i] = 0;
        p[i] = r[i];
    }
    ro = alpha = omega = 1.0;

    complex<double> a1 = 0.0, a2 = 0.0, a3;

    bool not_end = true;
    double discr = 2.0 * sqrt(dot_prod(r, r).real());
    size_t iter;
    for(iter = 0; iter < max_iter && not_end; iter++)
    {
        discr = sqrt(dot_prod(r, r).real());

        if(iter%10 == 0)
        {
            printf("BiCGStabVC Residual:\t%5lu\t%.3e\r", (unsigned long)iter, discr / rp_norm);
            fflush(stdout);
        }

        if(discr / rp_norm > eps)
        {
            mul_matrix(p, v);
            a1 = dot_prod_nocj(r, r2);
            a2 = dot_prod_nocj(v, r2);
            alpha = a1 / a2;
            for(size_t i = 0; i < n; i++)
                s[i] = r[i] - alpha * v[i];

            mul_matrix(s, t);
            a2 = dot_prod_nocj(s, t);
            a3 = dot_prod_nocj(t, t);
            omega = a2 / a3;

            for(size_t i = 0; i < n; i++)
            {
                x0[i] += alpha * p[i] + omega * s[i];
                r[i] = s[i] - omega * t[i];
            }

            a2 = dot_prod_nocj(r, r2);
            beta = a2 / a1 * alpha / omega;

            for(size_t i = 0; i < n; i++)
                p[i] = r[i] + beta * (p[i] - omega * v[i]);
        }
        else
            not_end = false;
    }

    mul_matrix(x0, r);
    for(size_t i = 0; i < n; i++)
        r[i] = rp[i] - r[i];
    discr = sqrt(dot_prod(r, r).real());
    printf("BiCGStabVC Residual:\t%5lu\t%.3e\n", (unsigned long)iter, discr / rp_norm);

    if(iter >= max_iter)
        printf("Soulution can`t found, iteration limit exceeded!\n");

    for(size_t i = 0; i < n; i++)
        solution[i] = x0[i];
    delete [] x0;
}
