#include "BiCGComplex_VC.h"
#if defined _MSC_VER && !defined _CRT_SECURE_NO_WARNINGS
#define _CRT_SECURE_NO_WARNINGS
#endif
#include <cstdio>
#include <cmath>

BiCGComplex_VC::BiCGComplex_VC()
{
    r = z = s = p = t = t1 = NULL;
}

BiCGComplex_VC::~BiCGComplex_VC()
{
    if(r) delete [] r;
    if(z) delete [] z;
    if(s) delete [] s;
    if(p) delete [] p;
    if(t) delete [] t;
    if(t1) delete [] t1;
}

void BiCGComplex_VC::init(size_t * gi_s, size_t * gj_s, complex<double> * di_s,
                          complex<double> * gg_s, size_t n_s)
{
    gi = gi_s;
    gj = gj_s;
    di = di_s;
    gg = gg_s;
    n = n_s;

    if(r) delete [] r;
    if(z) delete [] z;
    if(s) delete [] s;
    if(p) delete [] p;
    if(t) delete [] t;
    if(t1) delete [] t1;

    r = new complex<double> [n];
    z = new complex<double> [n];
    s = new complex<double> [n];
    p = new complex<double> [n];
    t = new complex<double> [n];
    t1 = new complex<double> [n];
}

complex<double> BiCGComplex_VC::dot_prod(const complex<double> * a, const complex<double> * b) const
{
    complex<double> res;
    for(size_t i = 0; i < n; i++)
    {
        res += conj(a[i]) * b[i];
    }
    return res;
}

complex<double> BiCGComplex_VC::dot_prod_nocj(const complex<double> * a, const complex<double> * b) const
{
    complex<double> res;
    for(size_t i = 0; i < n; i++)
    {
        res += a[i] * b[i];
    }
    return res;
}

void BiCGComplex_VC::mul_matrix(const complex<double> * f, complex<double> * x) const
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

void BiCGComplex_VC::solve(complex<double> * solution, complex<double> * rp, double gamma)
{
    double eps = gamma;
    size_t max_iter = /*(size_t) 100 * sqrt(n)*/ 15000;

    complex<double> dp1, dp2, alpha, beta;
    double rp_norm = sqrt(dot_prod(rp, rp).real());
    x0 = new complex<double> [n];

    mul_matrix(solution, t);
    for(size_t i = 0; i < n; i++)
    {
        x0[i] = solution[i];
        r[i] = rp[i] - t[i];
        p[i] = z[i] = s[i] = r[i];
    }

    bool not_end = true;
    double discr;
    size_t iter;

    dp1 = dot_prod_nocj(p, r);

    for(iter = 0; iter < max_iter && not_end; iter++)
    {
        discr = sqrt(dot_prod(r, r).real());

        if(iter%10 == 0)
        {
            printf("BiCGVC Residual:\t%5lu\t%.3e\r", (unsigned long)iter, discr / rp_norm);
            fflush(stdout);
        }

        if(discr / rp_norm > eps)
        {
            mul_matrix(z, t);

            alpha = dp1 / dot_prod_nocj(s, t);

            mul_matrix(s, t1);

            for(size_t i = 0; i < n; i++)
            {
                x0[i] += alpha * z[i];
                r[i] -= alpha * t[i];
                p[i] -= alpha * t1[i];
            }

            dp2 = dot_prod_nocj(p, r);
            beta = dp2 / dp1;
            dp1 = dp2;

            for(size_t i = 0; i < n; i++)
            {
                z[i] = r[i] + beta * z[i];
                s[i] = p[i] + beta * s[i];
            }
        }
        else
            not_end = false;
    }

    mul_matrix(x0, r);
    for(size_t i = 0; i < n; i++)
        r[i] = rp[i] - r[i];
    discr = sqrt(dot_prod(r, r).real());
    printf("BiCGVC Residual:\t%5lu\t%.3e\t%.3e\n", (unsigned long)iter, discr / rp_norm, gamma);

    if(iter >= max_iter)
        printf("Soulution can`t found, iteration limit exceeded!\n");

    for(size_t i = 0; i < n; i++)
        solution[i] = x0[i];
    delete [] x0;
}
