#if defined _MSC_VER && !defined _CRT_SECURE_NO_WARNINGS
#define _CRT_SECURE_NO_WARNINGS
#endif
#include "GMRES_Complex_LLT.h"
#include <cstdio>
#include <cmath>

void GMRES_Complex_LLT::init(const size_t * gi_s, const size_t * gj_s, const complex<double> * di_s,
                             const complex<double> * gg_s, size_t n_s)
{
    gi = gi_s;
    gj = gj_s;
    di = di_s;
    gg = gg_s;
    n = n_s;

    delete [] r;
    delete [] w;
    delete [] t;
    delete [] L_di;
    delete [] L_gg;
    delete [] d;
    if(H && m)
        for(size_t i = 0; i <= m; i++)
            delete [] H[i];
    delete [] H;
    if(VT && m)
        for(size_t i = 0; i < m; i++)
            delete [] VT[i];
    delete [] VT;

    r = new complex<double> [n];
    w = new complex<double> [n];
    t = new complex<double> [n];
    L_di = new complex<double> [n];
    L_gg = new complex<double> [gi[n]];

    for(size_t i = 0; i < n; i++)
    {
        L_di[i] = di[i];
        //x0[i] = 0.0;
    }

    for(size_t i = 0 ; i < gi[n] ; i++)
        L_gg[i] = gg[i];

    make_LLT_decomposition();

    // Глубина метода
    m = m_curr = 5;
//    FILE * fp = fopen("gmres_m.txt", "r");
//    if(fp)
//    {
//        unsigned int _m;
//        fscanf(fp, "%u", & _m);
//        m = m_curr = _m;
//        fclose(fp);
//    }
    char * env_m = getenv("GMRES_M");
    if(env_m)
    {
        int temp_m = atoi(env_m);
        if(temp_m >= 1)
            m = m_curr = (size_t)temp_m;
    }

    VT = new complex<double> * [m];
    for(size_t i = 0; i < m; i++)
        VT[i] = new complex<double> [n];
    H = new complex<double> * [m + 1];
    for(size_t i = 0; i < m + 1; i++)
    {
        H[i] = new complex<double> [m];
        for(size_t j = 0; j < m; j++)
            H[i][j] = 0;
    }
    d = new complex<double> [m + 1];
}

void GMRES_Complex_LLT::make_LLT_decomposition()
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
                        sum_l += L_gg[m] * L_gg[j];
                        j_s++;
                    }
                }
            }
            L_gg[i] = (L_gg[i] -  sum_l) / L_di[gj[i]];

            sum_d += L_gg[i] * L_gg[i];
        }
        L_di[k] = sqrt(L_di[k] - sum_d);
    }
}

complex<double> GMRES_Complex_LLT::dot_prod(const complex<double> * a, const complex<double> * b) const
{
    complex<double> d_p = 0.0;
    for(size_t i = 0; i < n; i++)
        d_p += conj(a[i]) * b[i];
    return d_p;
}

double GMRES_Complex_LLT::dot_prod_self(const complex<double> * a) const
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

void GMRES_Complex_LLT::mul_matrix(const complex<double> * f, complex<double> * x) const
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

void GMRES_Complex_LLT::solve_L(const complex<double> * f, complex<double> * x) const
{
    for(size_t k = 1, k1 = 0; k <= n; k++, k1++)
    {
        complex<double> sum = 0.0;

        for(size_t i = gi[k1]; i < gi[k]; i++)
            sum += L_gg[i] * x[gj[i]];

        x[k1] = (f[k1] - sum) / L_di[k1];
    }
}

void GMRES_Complex_LLT::solve_LT(complex<double> * f, complex<double> * x) const
{
    for(size_t k = n, k1 = n-1; k > 0; k--, k1--)
    {
        x[k1] = f[k1] / L_di[k1];
        complex<double> v_el = x[k1];

        for(size_t i = gi[k1]; i < gi[k]; i++)
            f[gj[i]] -= L_gg[i] * v_el;
    }
}

void GMRES_Complex_LLT::solve_LTAL(const complex<double> * f, complex<double> * x, complex<double> * tmp) const
{
    copy_vec(f, x);
    solve_LT(x, x);
    mul_matrix(x, tmp);
    solve_L(tmp, x);
}

void GMRES_Complex_LLT::mul_LT(const complex<double> * f, complex<double> * x) const
{
    for(size_t i = 0; i < n; i++)
    {
        complex<double> v_el = f[i];
        x[i] = L_di[i] * v_el;
        for(size_t k = gi[i], k1 = gi[i+1]; k < k1; k++)
            x[gj[k]] += L_gg[k] * v_el;
    }
}

void GMRES_Complex_LLT::copy_vec(const complex<double> * f, complex<double> * x) const
{
    for(size_t i = 0; i < n; i++)
        x[i] = f[i];
}

bool GMRES_Complex_LLT::is_fpu_error(double x) const
{
    double y = x - x;
    return x != x || y != y;
}

void GMRES_Complex_LLT::solve(complex<double> * solution, const complex<double> * rp_s,
                              double eps, size_t max_iter)
{
    rp = rp_s;

    x0 = new complex<double> [n];
    mul_LT(solution, x0);

    mul_matrix(solution, r);
    for(size_t i = 0; i < n ; i++)
        t[i] = rp[i] - r[i];
    solve_L(t, r);

    solve_L(rp, w); // Правая часть у нас предобусловленная
    double discr = 0.0, rp_norm = sqrt(dot_prod_self(w));
    double residual = 0.0, residual_prev;

    if(is_fpu_error(rp_norm))
    {
        fprintf(stderr, "Error: FPU error detected in right part!\n");
        solve_LT(x0, solution);
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
            solve_LT(x0, solution);
            delete [] x0;
            return;
        }

        residual_prev = residual;
        residual = discr / rp_norm;
        if(iter%10 == 0)
        {
            printf("GMRES_Complex_LLT Residual:\t%5lu\t%.3e\r", (unsigned long)iter, residual);
            fflush(stdout);
        }

        if(residual > eps && fabs(residual - residual_prev) / (residual) > 1e-15)
        {
            // V
            for(size_t i = 0; i < n; i++)
                VT[0][i] = r[i] / discr;

            // d
            d[0] = discr;
            for(size_t i = 1; i <= m_curr; i++)
                d[i] = 0.0;

            for(size_t mu = 0; mu < m_curr; mu++)
            {
                // w
                solve_LTAL(VT[mu], w, t);

                // H_lambda_mu
                for(size_t lambda = 0; lambda <= mu; lambda++)
                    H[lambda][mu] = dot_prod(VT[lambda], w);

                // tilde_v_mu+1
                for(size_t lambda = 0; lambda <= mu; lambda++)
                    for(size_t i = 0; i < n; i++)
                        w[i] -= VT[lambda][i] * H[lambda][mu];

                // H_mu+1_mu
                double vnorm = sqrt(dot_prod_self(w));
                H[mu + 1][mu] = vnorm;

                if(abs(vnorm) < 1e-15)
                {
                    if(mu > 1)  m_curr = mu;
                    else        m_curr = 1;
                    break;
                }

                // v_mu+1
                if(mu + 1 < m_curr)
                    for(size_t i = 0; i < n; i++)
                        VT[mu + 1][i] = w[i] / vnorm;
            }

            // H(i), d(i)
            for(size_t i = 0; i < m_curr; i++)
            {
                double t1 = H[i + 1][i].real();
                double t2 = abs(H[i][i]);
                double t3 = 1.0 / sqrt(t2 * t2 + t1 * t1);
                double s = t1 * t3;
                complex<double> c = H[i][i] * t3;

                for(size_t j = i; j < m_curr; j++)
                {
                    complex<double> H_ij = H[i][j];
                    H[i][j] = H_ij * conj(c) + H[i + 1][j] * s;
                    H[i + 1][j] = - H_ij * s + H[i + 1][j] * c;
                }
                complex<double> d_i = d[i];
                d[i] = d_i * conj(c) + d[i + 1] * s;
                d[i + 1] = - d_i * s + d[i + 1] * c;
            }

            // Hz = d
            for(size_t i = 0; i < m_curr; i++)
            {
                d[i] = d[i] / H[i][i];
                for(size_t j = i + 1; j < m_curr; j++)
                    H[i][j] /= H[i][i];
            }
            for(size_t i = m_curr - 1; i > 0; i--)
                for(size_t j = i - 1, ju = i; ju > 0; j--, ju--)
                    d[j] -= H[j][i] * d[i];

            // x = x + Vz
            for(size_t i = 0; i < n; i++)
                for(size_t j = 0; j < m_curr; j++)
                    x0[i] += VT[j][i] * d[j];

            // r
            copy_vec(x0, t);
            solve_LT(t, t);
            mul_matrix(t, r);
            for(size_t i = 0; i < n; i++)
                t[i] = rp[i] - r[i];
            solve_L(t, r);
        }
        else
            finished = true;
    }

    solve_LT(x0, solution);
    delete [] x0;

//    mul_matrix(solution, r);
//    for(size_t i = 0; i < n; i++)
//        r[i] = rp[i] - r[i];
//    discr = sqrt(dot_prod_self(r));
//    rp_norm = sqrt(dot_prod_self(rp));
    printf("GMRES_Complex_LLT Residual:\t%5lu\t%.3e\n", (unsigned long)iter - 1, discr / rp_norm);

    if(iter >= max_iter)
        printf("Soulution can`t found, iteration limit exceeded!\n");
    else if(residual > eps)
        printf("Soulution can`t found, stagnation detected!\n");
}

GMRES_Complex_LLT::GMRES_Complex_LLT()
{
    d = r = w = t = x0 = L_di = L_gg = NULL;
    H = VT = NULL;
    m = m_curr = 0;
}

GMRES_Complex_LLT::~GMRES_Complex_LLT()
{
    delete [] r;
    delete [] w;
    delete [] t;
    delete [] L_di;
    delete [] L_gg;
    delete [] d;
    if(H && m)
        for(size_t i = 0; i <= m; i++)
            delete [] H[i];
    delete [] H;
    if(VT && m)
        for(size_t i = 0; i < m; i++)
            delete [] VT[i];
    delete [] VT;
}
