#if defined _MSC_VER && !defined _CRT_SECURE_NO_WARNINGS
#define _CRT_SECURE_NO_WARNINGS
#endif
#include "V_cycle.h"
#include <cstdio>
#include <cmath>

void V_cycle::init_lvl1(size_t * gi_s, size_t * gj_s, complex<double> * di_s, complex<double> * gg_s, complex<double> * rp_s, size_t n_h_s)
{
    n_lvl1 = n_h_s;

    gi_lvl1 = gi_s;
    gj_lvl1 = gj_s;
    di_lvl1 = di_s;
    gg_lvl1 = gg_s;
    rp = rp_s;
    bcgm_lvl1.init(gi_s, gj_s, di_s, gg_s, n_lvl1);
}

void V_cycle::init_lvl2(size_t * gi_s, size_t * gj_s, complex<double> * di_s, complex<double> * gg_s, size_t n_2h_s)
{
    n_lvl2 = n_2h_s;
    bcgm_lvl2.init(gi_s, gj_s, di_s, gg_s, n_lvl2);
}

void V_cycle::init_operators(size_t * gi_R_s, size_t * gj_R_s, double * gg_R_s)
{
    gi_R = gi_R_s;
    gj_R = gj_R_s;
    gg_R = gg_R_s;
}

void V_cycle::get_grad_bounds(set<size_t> * grad_bounds_s)
{
    grad_bounds = grad_bounds_s;
}

void V_cycle::get_main_bounds(set<size_t> * main_bounds_s)
{
    main_bounds = main_bounds_s;
}


complex<double> V_cycle::dot_prod(complex<double> * a, complex<double> * b)
{
    complex<double> d_p = 0;
    for(size_t i = 0; i < n_lvl1; i++)
        d_p += conj(a[i]) * b[i];
    return d_p;
}

void V_cycle::mul_matrix(complex<double> * f, complex<double> *& x)
{

    for(size_t i = 0; i < n_lvl1; i++)
    {
        complex<double> v_el = f[i];
        x[i] = di_lvl1[i] * v_el;
        for(size_t k = gi_lvl1[i], k1 = gi_lvl1[i + 1]; k < k1; k++)
        {
            size_t j = gj_lvl1[k];
            x[i] += gg_lvl1[k] * f[j];
            x[j] += gg_lvl1[k] * v_el;
        }
    }
}

void V_cycle::calc_nev(complex<double> * x0, complex<double> * p)
{
    mul_matrix(x0, p);
    for(size_t i = 0; i < n_lvl1; i++)
        p[i] = rp[i] - p[i];

}

void V_cycle::solve(complex<double> * solution, double eps)
{
    size_t max_iter = 1000;

    complex<double> * x0 = new complex<double> [n_lvl1];
    complex<double> * d = new complex<double> [n_lvl2];
    complex<double> * p = new complex<double> [n_lvl1];
    complex<double> * r_n = new complex<double> [n_lvl1];
    complex<double> * e_2h = new complex<double> [n_lvl2];

    for(size_t i = 0 ; i < n_lvl1; i++)
        x0[i] = solution[i];

    for(size_t i = 0; i < n_lvl2; i++)
        e_2h[i] = 0;

    bool end = false;

    double rp_norm = sqrt(dot_prod(rp, rp).real());
    size_t iter;
    calc_nev(x0, p);
    for(iter = 0; iter < max_iter && !end; iter++)
    {
        for(size_t i = 0; i < n_lvl2; i++)
        {
            d[i] = 0;
            for(size_t j = gi_R[i]; j < gi_R[i + 1]; j++)
                d[i] += gg_R[j] * p[gj_R[j]];
        }
        for(set<size_t>::const_iterator i = grad_bounds->begin(); i != grad_bounds->end(); i++)
            d[*i] = 0;

        bcgm_lvl2.solve(e_2h, d, 0.05);

        for(size_t i = 0; i < n_lvl1; i++)
            r_n[i] = 0;

        for(size_t i = 0; i < n_lvl2; i++)
        {
            for(size_t j = gi_R[i]; j < gi_R[i + 1]; j++)
                r_n[gj_R[j]] += gg_R[j] * e_2h[i];
        }

        //Уточнение на градиентном пространсте
        for(size_t i = 0; i < n_lvl1; i++)
            x0[i] += r_n[i];

        calc_nev(x0, p);

        for(set<size_t>::const_iterator it = main_bounds->begin(); it != main_bounds->end(); it++)
            p[(*it)] = 0;

        bcgm_lvl1.solve(r_n, p, 0.1);

        //Уточненение на всём пространстве
        for(size_t i = 0; i < n_lvl1; i++)
            x0[i] += r_n[i];

        calc_nev(x0, p);
        double e_h = sqrt(dot_prod(p,p).real());

        if(e_h / rp_norm < eps)
            end = true;
        //if(iter%50 == 0)
            printf("%d\t%3e\n", (unsigned)iter, e_h / rp_norm);
    }

    for(size_t i = 0; i < n_lvl1; i++)
        solution[i] = x0[i];

    delete [] x0;
    delete [] p;
    delete [] d;
    delete [] r_n;
    delete [] e_2h;
}
