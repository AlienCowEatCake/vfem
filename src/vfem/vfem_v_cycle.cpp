#include "vfem.h"

// Проектирование на пространство ядра
void VFEM::to_kernel_space(const complex<double> * in, complex<double> * out) const
{
    size_t nodes_num = nodes.size();
    MAYBE_UNUSED(nodes_num);
    size_t edges_num = edges.size();
    MAYBE_UNUSED(edges_num);
#if BASIS_ORDER >= 2
    size_t faces_num = faces.size();
    MAYBE_UNUSED(faces_num);
#endif

    for(size_t i = 0; i < ker_dof_num; i++)
        out[i] = 0.0;

    for(set<edge>::iterator it = edges.begin(); it != edges.end(); ++it)
    {
        out[it->nodes[0]->num] -= in[it->num];
        out[it->nodes[1]->num] += in[it->num];
    }

#if BASIS_ORDER >= 2 || BASIS_TYPE == 2
    for(set<edge>::iterator it = edges.begin(); it != edges.end(); ++it)
        out[nodes_num + it->num] = in[edges_num + it->num];
#endif

#if BASIS_ORDER > 2 || (BASIS_TYPE == 2 && BASIS_ORDER == 2)
    for(set<face>::iterator it = faces.begin(); it != faces.end(); ++it)
        out[nodes_num + edges_num + it->num] = in[2 * edges_num + 2 * faces_num + it->num];
    for(set<edge>::iterator it = edges.begin(); it != edges.end(); ++it)
        out[nodes_num + edges_num + faces_num + it->num] = in[2 * edges_num + 3 * faces_num + it->num];
#endif
}

// Интерполяция на полное пространство
void VFEM::to_full_space(const complex<double> * in, complex<double> * out) const
{
    size_t nodes_num = nodes.size();
    MAYBE_UNUSED(nodes_num);
    size_t edges_num = edges.size();
    MAYBE_UNUSED(edges_num);
#if BASIS_ORDER >= 2
    size_t faces_num = faces.size();
    MAYBE_UNUSED(faces_num);
#endif

    for(size_t i = 0; i < dof_num; i++)
        out[i] = 0.0;

    for(set<edge>::iterator it = edges.begin(); it != edges.end(); ++it)
    {
        out[it->num] -= in[it->nodes[0]->num];
        out[it->num] += in[it->nodes[1]->num];
    }

#if BASIS_ORDER >= 2 || BASIS_TYPE == 2
    for(set<edge>::iterator it = edges.begin(); it != edges.end(); ++it)
        out[edges_num + it->num] = in[nodes_num + it->num];
#endif

#if BASIS_ORDER > 2 || (BASIS_TYPE == 2 && BASIS_ORDER == 2)
    for(set<face>::iterator it = faces.begin(); it != faces.end(); ++it)
        out[2 * edges_num + 2 * faces_num + it->num] = in[nodes_num + edges_num + it->num];
    for(set<edge>::iterator it = edges.begin(); it != edges.end(); ++it)
        out[2 * edges_num + 3 * faces_num + it->num] = in[nodes_num + edges_num + faces_num + it->num];
#endif
}

// Скалярное произведение
double VFEM::dot_prod_self(const complex<double> * a) const
{
    double d_p = 0.0;
    for(size_t i = 0; i < dof_num; i++)
    {
        double re = a[i].real();
        double im = a[i].imag();
        d_p += re * re + im * im;
    }
    return d_p;
}

// Умножение матрицы с полного пространства на вектор
void VFEM::mul_matrix(const complex<double> * f, complex<double> * x) const
{
    for(size_t i = 0; i < dof_num; i++)
    {
        complex<double> v_el = f[i];
        x[i] = slae.di[i] * v_el;
        for(size_t k = slae.ig[i], k1 = slae.ig[i + 1]; k < k1; k++)
        {
            size_t j = slae.jg[k];
            x[i] += slae.gg[k] * f[j];
            x[j] += slae.gg[k] * v_el;
        }
    }
}

// Подсчет невязки
void VFEM::calc_residual(const complex<double> * x0, complex<double> * p) const
{
    mul_matrix(x0, p);
    for(size_t i = 0; i < dof_num; i++)
        p[i] = slae.rp[i] - p[i];
}

// Запуск решения СЛАУ
void VFEM::solve()
{
    extern double SLAE_MAIN_EPSILON;
    //slae.solve(SLAE_MAIN_EPSILON);
    //return;

    size_t max_iter = 1000;
    double gamma0       = 0.1;
    double gamma_full   = 0.2;
    double gamma_ker    = 0.3;

    slae.inline_init();
    ker_slae.inline_init();

    // Уточнение начального приближения на полном пространстве
    slae.inline_solve(slae.x, slae.rp, gamma0);

    double rp_norm = sqrt(dot_prod_self(slae.rp));

    // Вектор невязки
    complex<double> * r = new complex<double> [dof_num];
    calc_residual(slae.x, r);

    complex<double> * g = new complex<double> [ker_dof_num];
    complex<double> * y = new complex<double> [dof_num];

    for(size_t iter = 0; iter < max_iter; iter++)
    {
        // g = Pr
        to_kernel_space(r, g);

        // Правим краевые
        for(set<size_t>::iterator it = ker_dof_first.begin(); it != ker_dof_first.end(); ++it)
            g[*it] = 0.0;

        // z = (PAPt)^-1 g или z = solve1(PAPt, g)
        for(size_t i = 0; i < ker_dof_num; i++) ker_slae.x[i] = 0.0;
        ker_slae.inline_solve(ker_slae.x, g, gamma_ker);

        // Уточнение на градиентном пространсте
        // x = x + Ptz
        to_full_space(ker_slae.x, y);
        for(size_t i = 0; i < dof_num; i++)
            slae.x[i] += y[i];

        // r = b - Ax
        calc_residual(slae.x, r);
        printf("V-Cycle[K] Residual:\t%5lu\t%.3e\n", (unsigned long)iter, sqrt(dot_prod_self(r)) / rp_norm);

        // Правим краевые
#if defined VFEM_USE_NONHOMOGENEOUS_FIRST
        for(map<size_t, size_t>::iterator it = global_to_local.begin(); it != global_to_local.end(); ++it)
            r[it->first] = 0.0;
#else
        for(set<size_t>::iterator it = dof_first.begin(); it != dof_first.end(); ++it)
            r[(*it)] = 0.0;
#endif

        // y = solve2(A, r)
        for(size_t i = 0; i < dof_num; i++) y[i] = 0.0;
        slae.inline_solve(y, r, gamma_full);

        //Уточненение на всём пространстве
        // x = x + y
        for(size_t i = 0; i < dof_num; i++)
            slae.x[i] += y[i];

        // r = b - Ax
        calc_residual(slae.x, r);
        double res = sqrt(dot_prod_self(r)) / rp_norm;
        printf("V-Cycle[F] Residual:\t%5lu\t%.3e\n", (unsigned long)iter, res);

        if(res < SLAE_MAIN_EPSILON) break;
    }

    delete [] r;
    delete [] g;
    delete [] y;
}
