#if defined(_MSC_VER) && !defined(_CRT_SECURE_NO_WARNINGS)
#define _CRT_SECURE_NO_WARNINGS
#endif
#include "COCR_LDLT.h"
#include <cstdio>
#include <cmath>

void COCR_LDLT::init(const size_t * gi_s, const size_t * gj_s, const complex<double> * di_s,
                     const complex<double> * gg_s, size_t n_s)
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
    delete [] w;
    delete [] a;
    delete [] L_di;
    delete [] L_gg;

    r = new complex<double> [n];
    z = new complex<double> [n];
    p = new complex<double> [n];
    s = new complex<double> [n];
    w = new complex<double> [n];
    a = new complex<double> [n];

    L_di = new complex<double> [n];
    L_gg = new complex<double> [m];

    for(size_t i = 0; i < n; i++)
    {
        L_di[i] = di[i];
        //x0[i] = 0.0;
    }

    for(size_t i = 0 ; i < m ; i++)
        L_gg[i] = gg[i];

    make_LDLT_decomposition();
}

void COCR_LDLT::make_LDLT_decomposition()
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
                        //sum_l += L_gg[m] * L_gg[j] * L_di[gj[m]];
                        sum_l += L_gg[m] * L_gg[j] / L_di[gj[m]]; // Warning: Инвертировано!
                        j_s++;
                    }
                }
            }
            //L_gg[i] = (L_gg[i] -  sum_l) / L_di[gj[i]];
            L_gg[i] = (L_gg[i] -  sum_l) * L_di[gj[i]]; // Warning: Инвертировано!

            //sum_d += L_gg[i] * L_gg[i] * L_di[gj[i]];
            sum_d += L_gg[i] * L_gg[i] / L_di[gj[i]]; // Warning: Инвертировано!
        }
        //L_di[k] = L_di[k] - sum_d;
        L_di[k] = 1.0 / (L_di[k] - sum_d); // Warning: Инвертировано!
    }
}

complex<double> COCR_LDLT::dot_prod_nocj(const complex<double> * a, const complex<double> * b) const
{
    complex<double> d_p = 0.0;
    for(size_t i = 0; i < n; i++)
        d_p += a[i] * b[i];
    return d_p;
}

double COCR_LDLT::dot_prod_self(const complex<double> * a) const
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

void COCR_LDLT::mul_matrix(const complex<double> * f, complex<double> * x) const
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

void COCR_LDLT::solve_L(const complex<double> * f, complex<double> * x) const
{
    for(size_t k = 1, k1 = 0; k <= n; k++, k1++)
    {
        complex<double> sum = 0.0;

        for(size_t i = gi[k1]; i < gi[k]; i++)
            sum += L_gg[i] * x[gj[i]];

        x[k1] = (f[k1] - sum);
    }
}

void COCR_LDLT::solve_LT(complex<double> * f, complex<double> * x) const
{
    for(size_t k = n, k1 = n-1; k > 0; k--, k1--)
    {
        x[k1] = f[k1];
        complex<double> v_el = x[k1];

        for(size_t i = gi[k1]; i < gi[k]; i++)
            f[gj[i]] -= L_gg[i] * v_el;
    }
}

void COCR_LDLT::solve_LDLT(const complex<double> * f, complex<double> * x) const
{
    solve_L(f, x);
    for(size_t i = 0; i < n; i++)
        //x[i] /= L_di[i];
        x[i] *= L_di[i]; // Warning: Инвертировано!
    solve_LT(x, x);
}

bool COCR_LDLT::is_fpu_error(double x) const
{
    double y = x - x;
    return x != x || y != y;
}

void COCR_LDLT::solve(complex<double> * solution, const complex<double> * rp_s,
                      double eps, size_t max_iter)
{
    rp = rp_s;

    x0 = new complex<double> [n];
    for(size_t i = 0; i < n; i++)
        x0[i] = solution[i];

    mul_matrix(x0, r);

    for(size_t i = 0; i < n ; i++)
        r[i] = rp[i] - r[i];

    solve_LDLT(r, s);
    mul_matrix(s, z);
    solve_LDLT(z, w);
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
            printf("COCR_LDLT Residual:\t%5lu\t%.3e\r", (unsigned long)iter, residual);
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
            }
            solve_LDLT(z, w);
        }
        else
            finished = true;
    }

//    mul_matrix(x0, r);
//    for(size_t i = 0; i < n; i++)
//        r[i] = rp[i] - r[i];
//    discr = sqrt(dot_prod_self(r));
    printf("COCR_LDLT Residual:\t%5lu\t%.3e\n", (unsigned long)iter - 1, discr / rp_norm);

    if(iter >= max_iter)
        printf("Soulution can`t found, iteration limit exceeded!\n");
    else if(residual > eps)
        printf("Soulution can`t found, stagnation detected!\n");

    for(size_t i = 0; i < n; i++)
        solution[i] = x0[i];
    delete [] x0;
}

COCR_LDLT::COCR_LDLT()
{
    r = x0 = z = p = s = L_di = L_gg = w = a = NULL;
}

COCR_LDLT::~COCR_LDLT()
{
    delete [] r;
    delete [] z;
    delete [] p;
    delete [] s;
    delete [] w;
    delete [] a;
    delete [] L_di;
    delete [] L_gg;
}
