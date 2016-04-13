#if defined(_MSC_VER) && !defined(_CRT_SECURE_NO_WARNINGS)
#define _CRT_SECURE_NO_WARNINGS
#endif
#include "BiCGStab_Complex_Smooth.h"
#include <cstdio>
#include <cmath>

BiCGStab_Complex_Smooth::BiCGStab_Complex_Smooth()
{
    r = v = s = p = t = r2 = xs = rs = NULL;
}

BiCGStab_Complex_Smooth::~BiCGStab_Complex_Smooth()
{
    delete [] r;
    delete [] v;
    delete [] s;
    delete [] p;
    delete [] t;
    delete [] r2;
    delete [] xs;
    delete [] rs;
}

void BiCGStab_Complex_Smooth::init(const size_t * gi_s, const size_t * gj_s, const complex<double> * di_s,
                              const complex<double> * gg_s, size_t n_s)
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
    delete [] xs;
    delete [] rs;

    r = new complex<double> [n];
    v = new complex<double> [n];
    s = new complex<double> [n];
    p = new complex<double> [n];
    t = new complex<double> [n];
    r2 = new complex<double> [n];
    xs = new complex<double> [n];
    rs = new complex<double> [n];
}

complex<double> BiCGStab_Complex_Smooth::dot_prod(const complex<double> * a, const complex<double> * b) const
{
    complex<double> res;
    for(size_t i = 0; i < n; i++)
    {
        res += conj(a[i]) * b[i];
    }
    return res;
}

complex<double> BiCGStab_Complex_Smooth::dot_prod_nocj(const complex<double> * a, const complex<double> * b) const
{
    complex<double> res;
    for(size_t i = 0; i < n; i++)
    {
        res += a[i] * b[i];
    }
    return res;
}

double BiCGStab_Complex_Smooth::dot_prod_self(const complex<double> * a) const
{
    double res = 0.0;
    for(size_t i = 0; i < n; i++)
    {
        double re = a[i].real();
        double im = a[i].imag();
        res += re * re + im * im;
    }
    return res;
}

double BiCGStab_Complex_Smooth::dot_prod_real(const complex<double> * a, const complex<double> * b) const
{
    double d_p = 0.0;
    for(size_t i = 0; i < n; i++)
        d_p += a[i].real() * b[i].real() + a[i].imag() * b[i].imag();
    return d_p;
}

void BiCGStab_Complex_Smooth::mul_matrix(const complex<double> * f, complex<double> * x) const
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

void BiCGStab_Complex_Smooth::solve(complex<double> * solution, const complex<double> * rp,
                                    double gamma, size_t max_iter)
{
    double eps = gamma * gamma;

    complex<double> omega, alpha, ro, beta;
    double rp_norm = dot_prod_self(rp);
    x0 = solution;

    mul_matrix(solution, t);
    for(size_t i = 0; i < n; i++)
    {
        xs[i] = solution[i];
        rs[i] = p[i] = r[i] = rp[i] - t[i];
        r2[i] = 1.0;
        v[i] = 0;
    }
    ro = alpha = omega = 1.0;

    complex<double> a1 = 0.0, a2 = 0.0, a3;

    bool not_end = true;
    double discr = 0.0;
    size_t iter;
    for(iter = 0; iter < max_iter && not_end; iter++)
    {
        discr = dot_prod_self(rs);

        if(iter%10 == 0)
        {
            printf("BiCGStab_Complex_Smooth Residual:\t%5lu\t%.3e\r", (unsigned long)iter, sqrt(discr / rp_norm));
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
                s[i] = r[i] - rs[i]; // r - s, сглаживатель
            }

            double eta = - dot_prod_real(rs, s) / dot_prod_self(s);
            if(eta < 0.0) eta = 0.0;
            else if(eta > 1.0) eta = 1.0;
            for(size_t i = 0; i < n ; i++)
            {
                double eta1 = 1.0 - eta;
                xs[i] = eta1 * xs[i] + eta * x0[i];
                rs[i] = eta1 * rs[i] + eta * r[i];
            }

            a2 = dot_prod_nocj(r, r2);
            beta = a2 / a1 * alpha / omega;

            for(size_t i = 0; i < n; i++)
                p[i] = r[i] + beta * (p[i] - omega * v[i]);
        }
        else
            not_end = false;
    }

//    mul_matrix(xs, r);
//    for(size_t i = 0; i < n; i++)
//        r[i] = rp[i] - r[i];
//    discr = dot_prod_self(r);
    printf("BiCGStab_Complex_Smooth Residual:\t%5lu\t%.3e\n", (unsigned long)iter - 1, sqrt(discr / rp_norm));

    if(iter >= max_iter)
        printf("Soulution can`t found, iteration limit exceeded!\n");

    for(size_t i = 0; i < n; i++)
        solution[i] = xs[i];
}
