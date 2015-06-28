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
    gi_lvl2 = gi_s;
    gj_lvl2 = gj_s;
    di_lvl2 = di_s;
    gg_lvl2 = gg_s;
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


complex<double> V_cycle::dot_prod(const complex<double> * a, const complex<double> * b) const
{
    complex<double> d_p = 0;
    for(size_t i = 0; i < n_lvl1; i++)
        d_p += conj(a[i]) * b[i];
    return d_p;
}

double V_cycle::dot_prod_self(const complex<double> * a) const
{
    double d_p = 0.0;
    for(size_t i = 0; i < n_lvl1; i++)
    {
        double re = a[i].real();
        double im = a[i].imag();
        d_p += re * re + im * im;
    }
    return d_p;
}

void V_cycle::mul_matrix(const complex<double> * f, complex<double> * x) const
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

void V_cycle::calc_residual(const complex<double> * x0, complex<double> * p) const
{
    mul_matrix(x0, p);
    for(size_t i = 0; i < n_lvl1; i++)
        p[i] = rp[i] - p[i];
}

void V_cycle::solve(complex<double> * solution, double eps)
{
    bool debug_print = false;
    if(debug_print)
    {
        size_t n = n_lvl2;
        size_t m = n_lvl1;
        double ** A = new double * [n];
        for(size_t i = 0; i < n; i++)
        {
            A[i] = new double [m];
            memset(A[i], 0, sizeof(double) * (m));
        }
        for(size_t i = 0; i < n; i++)
        {
            for(size_t k = gi_R[i]; k < gi_R[i + 1]; k++)
            {
                A[i][gj_R[k]] = gg_R[k];
            }
        }

        ofstream ofs("R.txt");
        for(size_t i = 0; i < n; i++)
        {
            for(size_t j = 0; j < m - 1; j++)
                ofs << A[i][j] << "\t";
            ofs << A[i][m - 1] << endl;
        }
        ofs.flush();
        ofs.close();
    }
    if(debug_print)
    {
        size_t n = n_lvl1;
        complex<double> ** A = new complex<double> * [n];
        for(size_t i = 0; i < n; i++)
        {
            A[i] = new complex<double> [n];
            memset(A[i], 0, sizeof(complex<double>) * (n));
        }
        for(size_t i = 0; i < n; i++)
        {
            for(size_t k = gi_lvl1[i]; k < gi_lvl1[i + 1]; k++)
            {
                A[i][gj_lvl1[k]] = gg_lvl1[k];
                A[gj_lvl1[k]][i] = gg_lvl1[k];
            }
            A[i][i] = di_lvl1[i];
        }

        ofstream ofs("Re_lvl1.txt");
        ofs.precision(17);
        ofs << scientific;
        for(size_t i = 0; i < n; i++)
        {
            for(size_t j = 0; j < n - 1; j++)
                ofs << A[i][j].real() << "\t";
            ofs << A[i][n - 1].real() << endl;
        }
        ofs.flush();
        ofs.close();

        ofs.open("Im_lvl1.txt");
        ofs.precision(17);
        ofs << scientific;
        for(size_t i = 0; i < n; i++)
        {
            for(size_t j = 0; j < n - 1; j++)
                ofs << A[i][j].imag() << "\t";
            ofs << A[i][n - 1].imag() << endl;
        }
        ofs.flush();
        ofs.close();
    }
    if(debug_print)
    {
        size_t n = n_lvl2;
        complex<double> ** A = new complex<double> * [n];
        for(size_t i = 0; i < n; i++)
        {
            A[i] = new complex<double> [n];
            memset(A[i], 0, sizeof(complex<double>) * (n));
        }
        for(size_t i = 0; i < n; i++)
        {
            for(size_t k = gi_lvl2[i]; k < gi_lvl2[i + 1]; k++)
            {
                A[i][gj_lvl2[k]] = gg_lvl2[k];
                A[gj_lvl2[k]][i] = gg_lvl2[k];
            }
            A[i][i] = di_lvl2[i];
        }

        ofstream ofs("Re_lvl2.txt");
        ofs.precision(17);
        ofs << scientific;
        for(size_t i = 0; i < n; i++)
        {
            for(size_t j = 0; j < n - 1; j++)
                ofs << A[i][j].real() << "\t";
            ofs << A[i][n - 1].real() << endl;
        }
        ofs.flush();
        ofs.close();

        ofs.open("Im_lvl2.txt");
        ofs.precision(17);
        ofs << scientific;
        for(size_t i = 0; i < n; i++)
        {
            for(size_t j = 0; j < n - 1; j++)
                ofs << A[i][j].imag() << "\t";
            ofs << A[i][n - 1].imag() << endl;
        }
        ofs.flush();
        ofs.close();
    }

    //bcgm_lvl1.solve(solution, rp, 0.01);
/*
    size_t max_iter = 1000;

    complex<double> * x0 = new complex<double> [n_lvl1];
    complex<double> * v_g = new complex<double> [n_lvl2];
    complex<double> * v_r0 = new complex<double> [n_lvl1];
    complex<double> * v_Py = new complex<double> [n_lvl1];
    complex<double> * v_y = new complex<double> [n_lvl2];

    for(size_t i = 0 ; i < n_lvl1; i++)
        x0[i] = solution[i];

    for(size_t i = 0; i < n_lvl2; i++)
        v_y[i] = 0;

    bool end = false;

    double rp_norm = sqrt(dot_prod(rp, rp).real());
    size_t iter;
    calc_residual(x0, v_r0);
    for(iter = 0; iter < max_iter && !end; iter++)
    {
        for(size_t i = 0; i < n_lvl2; i++)
        {
            v_g[i] = 0;
            for(size_t j = gi_R[i]; j < gi_R[i + 1]; j++)
                v_g[i] += gg_R[j] * v_r0[gj_R[j]];
        }
        for(set<size_t>::const_iterator i = grad_bounds->begin(); i != grad_bounds->end(); i++)
            v_g[*i] = 0;

        printf("Level 2:\n");
        bcgm_lvl2.solve(v_y, v_g, 0.1);

        for(size_t i = 0; i < n_lvl1; i++)
            v_Py[i] = 0;

        for(size_t i = 0; i < n_lvl2; i++)
        {
            for(size_t j = gi_R[i]; j < gi_R[i + 1]; j++)
                v_Py[gj_R[j]] += gg_R[j] * v_y[i];
        }

        //Уточнение на градиентном пространсте
        for(size_t i = 0; i < n_lvl1; i++)
            x0[i] += v_Py[i];

        calc_residual(x0, v_r0);
        printf("%d\t%3e\n", (unsigned)iter, sqrt(dot_prod(v_r0,v_r0).real()) / rp_norm);

        for(set<size_t>::const_iterator it = main_bounds->begin(); it != main_bounds->end(); it++)
            v_r0[(*it)] = 0;

        printf("Level 1:\n");
        bcgm_lvl1.solve(v_Py, v_r0, 0.9);

        //Уточненение на всём пространстве
        for(size_t i = 0; i < n_lvl1; i++)
            x0[i] += v_Py[i];

        calc_residual(x0, v_r0);
        double e_h = sqrt(dot_prod(v_r0,v_r0).real());

        if(e_h / rp_norm < eps)
            end = true;
        //if(iter%50 == 0)
            printf("%d\t%3e\n", (unsigned)iter, e_h / rp_norm);
    }

    for(size_t i = 0; i < n_lvl1; i++)
        solution[i] = x0[i];

    delete [] x0;
    delete [] v_r0;
    delete [] v_g;
    delete [] v_Py;
    delete [] v_y;
*/
    size_t max_iter = 1000;
    double gamma0 = 0.05;
    double gamma1 = 0.09;
    double gamma2 = 0.01;

    // Уточнение начального приближения на полном пространстве
    bcgm_lvl1.solve(solution, rp, gamma0);

    double rp_norm = sqrt(dot_prod_self(rp));

    // Вектор невязки
    complex<double> * r = new complex<double> [n_lvl1];
    calc_residual(solution, r);

    complex<double> * g = new complex<double> [n_lvl2];
    complex<double> * z = new complex<double> [n_lvl2];
    complex<double> * y = new complex<double> [n_lvl1];

    for(size_t iter = 0; iter < max_iter; iter++)
    {
        // g = Pr
        for(size_t i = 0; i < n_lvl2; i++)
            g[i] = 0.0;
        for(size_t i = 0; i < n_lvl2; i++)
        {
            for(size_t j = gi_R[i]; j < gi_R[i + 1]; j++)
                g[i] += gg_R[j] * r[gj_R[j]];
        }
        // Правим краевые
        for(set<size_t>::iterator it = grad_bounds->begin(); it != grad_bounds->end(); ++it)
            g[*it] = 0.0;

        // z = (PAPt)^-1 g или z = solve1(PAPt, g)
        for(size_t i = 0; i < n_lvl2; i++) z[i] = 0.0;
        bcgm_lvl2.solve(z, g, gamma2);

        // Уточнение на градиентном пространсте
        // x = x + Ptz
        for(size_t i = 0; i < n_lvl2; i++)
            for(size_t j = gi_R[i]; j < gi_R[i + 1]; j++)
                solution[gj_R[j]] += gg_R[j] * z[i];

        // r = b - Ax
        calc_residual(solution, r);
        printf("[%u] Kersel space completed, residual = %3e\n", (unsigned)iter, sqrt(dot_prod_self(r)) / rp_norm);

        // Правим краевые
        for(set<size_t>::iterator it = main_bounds->begin(); it != main_bounds->end(); ++it)
            r[(*it)] = 0.0;

        // y = solve2(A, r)
        for(size_t i = 0; i < n_lvl1; i++) y[i] = 0.0;
        bcgm_lvl1.solve(y, r, gamma1);

        //Уточненение на всём пространстве
        // x = x + y
        for(size_t i = 0; i < n_lvl1; i++)
            solution[i] += y[i];

        // r = b - Ax
        calc_residual(solution, r);
        double res = sqrt(dot_prod_self(r)) / rp_norm;
        printf("[%u] Full space completed, residual = %3e\n", (unsigned)iter, res);

        if(res < eps) break;
    }

    delete [] r;
    delete [] g;
    delete [] z;
    delete [] y;
}
