#include "CGMComplex_VC.h"
#if defined _MSC_VER && !defined _CRT_SECURE_NO_WARNINGS
#define _CRT_SECURE_NO_WARNINGS
#endif
#include <cstdio>
#include <cmath>

CGMComplex_VC::CGMComplex_VC()
{
    r = w = s = p = t = u = NULL;
}

CGMComplex_VC::~CGMComplex_VC()
{
    delete [] r;
    delete [] w;
    delete [] s;
    delete [] p;
    delete [] t;
    delete [] u;
}

void CGMComplex_VC::init(size_t * gi_s, size_t * gj_s, complex<double> * di_s,
                         complex<double> * gg_s, size_t n_s)
{
    gi = gi_s;
    gj = gj_s;
    di = di_s;
    gg = gg_s;
    n = n_s;

    delete [] r;
    delete [] w;
    delete [] s;
    delete [] p;
    delete [] t;
    delete [] u;

    r = new complex<double> [n];
    w = new complex<double> [n];
    s = new complex<double> [n];
    p = new complex<double> [n];
    t = new complex<double> [n];
    u = new complex<double> [n];
}

complex<double> CGMComplex_VC::dot_prod(const complex<double> * a, const complex<double> * b) const
{
    complex<double> res;
    for(size_t i = 0; i < n; i++)
    {
        res += conj(a[i]) * b[i];
    }
    return res;
}

complex<double> CGMComplex_VC::dot_prod_nocj(const complex<double> * a, const complex<double> * b) const
{
    complex<double> res;
    for(size_t i = 0; i < n; i++)
    {
        res += a[i] * b[i];
    }
    return res;
}

void CGMComplex_VC::mul_matrix(const complex<double> * f, complex<double> * x) const
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

void CGMComplex_VC::solve(complex<double> * solution, complex<double> * rp, double gamma)
{
    double eps = gamma;
    size_t max_iter = /*(size_t)sqrt(n)*/ 15000;

    complex<double> alpha, beta, alpha1, alpha2;
    double rp_norm = sqrt(dot_prod(rp, rp).real());
    x0 = new complex<double> [n];

    mul_matrix(solution, t);
    for(size_t i = 0; i < n; i++)
    {
        x0[i] = solution[i];
        r[i] = rp[i] - t[i];
        w[i] = r[i] / di[i];
        p[i] = 0.0;
    }

    alpha1 = dot_prod_nocj(r, w);
    beta = 0.0;

    bool not_end = true;
    double discr;

    size_t iter;
    for(iter = 0; iter < max_iter && not_end; iter++)
    {
        discr = sqrt(dot_prod(r, r).real());

        if(iter%10 == 0)
        {
            printf("CCGMVC Residual:\t%5lu\t%.3e\r", (unsigned long)iter, discr / rp_norm);
            fflush(stdout);
        }

        if(discr / rp_norm > eps)
        {
            for(size_t i = 0; i < n; i++)
                p[i] = w[i] + beta * p[i];
            mul_matrix(p, u);
            alpha2 = dot_prod_nocj(u, p);
            alpha = alpha1 / alpha2;

            for(size_t i = 0; i < n; i++)
            {
                x0[i] += alpha * p[i];
                r[i] -= alpha * u[i];
                w[i] = r[i] / di[i];
            }

            alpha2 = alpha1;
            alpha1 = dot_prod_nocj(r, w);
            beta = alpha1 / alpha2;
        }
        else
            not_end = false;
    }

    mul_matrix(x0, r);
    for(size_t i = 0; i < n; i++)
        r[i] = rp[i] - r[i];
    discr = sqrt(dot_prod(r, r).real());
    printf("CCGMVC Residual:\t%5lu\t%.3e\n", (unsigned long)iter, discr / rp_norm);

    if(iter >= max_iter)
        printf("Soulution can`t found, iteration limit exceeded!\n");

    for(size_t i = 0; i < n; i++)
        solution[i] = x0[i];
    delete [] x0;
}
