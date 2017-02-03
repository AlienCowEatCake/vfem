#if defined(_MSC_VER) && !defined(_CRT_SECURE_NO_WARNINGS)
#define _CRT_SECURE_NO_WARNINGS
#endif
#include "GMRES_Complex.h"
#include <cstdio>
#include <cmath>
#include "../../../../../utils/fpu.h"

namespace fem_core { namespace solvers { namespace CSLR { namespace symmetric { namespace complex {

void GMRES_Complex::init(const std::size_t * gi, const std::size_t * gj, const std::complex<double> * di,
                         const std::complex<double> * gg, std::size_t n,
                         preconditioners::preconditioner_interface<std::complex<double>, std::size_t> * precond)
{
    m_gi = gi;
    m_gj = gj;
    m_di = di;
    m_gg = gg;
    m_n = n;
    m_precond = precond;

    delete [] m_r;
    delete [] m_w;
    delete [] m_t;
    delete [] m_d;
    if(m_H && m_m)
        for(std::size_t i = 0; i <= m_m; i++)
            delete [] m_H[i];
    delete [] m_H;
    if(m_VT && m_m)
        for(std::size_t i = 0; i < m_m; i++)
            delete [] m_VT[i];
    delete [] m_VT;

    m_r = new std::complex<double> [m_n];
    m_w = new std::complex<double> [m_n];
    m_t = new std::complex<double> [m_n];

    // Глубина метода
    m_m = m_m_curr = 5;
    char * env_m = getenv("GMRES_M");
    if(env_m)
    {
        int temp_m = atoi(env_m);
        if(temp_m >= 1)
            m_m = m_m_curr = static_cast<std::size_t>(temp_m);
    }

    m_VT = new std::complex<double> * [m_m];
    for(std::size_t i = 0; i < m_m; i++)
        m_VT[i] = new std::complex<double> [m_n];
    m_H = new std::complex<double> * [m_m + 1];
    for(std::size_t i = 0; i < m_m + 1; i++)
    {
        m_H[i] = new std::complex<double> [m_m];
        for(std::size_t j = 0; j < m_m; j++)
            m_H[i][j] = 0;
    }
    m_d = new std::complex<double> [m_m + 1];
}

std::complex<double> GMRES_Complex::dot_prod(const std::complex<double> * a, const std::complex<double> * b) const
{
    std::complex<double> d_p = 0.0;
    for(std::size_t i = 0; i < m_n; i++)
        d_p += conj(a[i]) * b[i];
    return d_p;
}

double GMRES_Complex::dot_prod_self(const std::complex<double> * a) const
{
    double d_p = 0.0;
    for(std::size_t i = 0; i < m_n; i++)
    {
        double re = a[i].real();
        double im = a[i].imag();
        d_p += re * re + im * im;
    }
    return d_p;
}

void GMRES_Complex::mul_matrix(const std::complex<double> * f, std::complex<double> * x) const
{
    for(std::size_t i = 0; i < m_n; i++)
    {
        std::complex<double> v_el = f[i];
        x[i] = m_di[i] * v_el;
        for(std::size_t k = m_gi[i], k1 = m_gi[i+1]; k < k1; k++)
        {
            std::size_t j = m_gj[k];
            x[i] += m_gg[k] * f[j];
            x[j] += m_gg[k] * v_el;
        }
    }
}

void GMRES_Complex::solve_QAS(const std::complex<double> * f, std::complex<double> * x, std::complex<double> * tmp) const
{
    copy_vec(f, x);
    m_precond->solve_Q(x, x);
    mul_matrix(x, tmp);
    m_precond->solve_S(tmp, x);
}

void GMRES_Complex::copy_vec(const std::complex<double> * f, std::complex<double> * x) const
{
    for(std::size_t i = 0; i < m_n; i++)
        x[i] = f[i];
}

void GMRES_Complex::solve(std::complex<double> * solution, const std::complex<double> * rp,
                          double eps, std::size_t max_iter)
{
    m_rp = rp;
    std::string precond_name_str = m_precond->get_name();
    const char * precond_name = precond_name_str.c_str();

    m_x0 = new std::complex<double> [m_n];
    m_precond->mul_Q(solution, m_x0);

    mul_matrix(solution, m_r);
    for(std::size_t i = 0; i < m_n ; i++)
        m_t[i] = m_rp[i] - m_r[i];
    m_precond->solve_S(m_t, m_r);

    m_precond->solve_S(m_rp, m_w); // Правая часть у нас предобусловленная
    double discr = 0.0, rp_norm = sqrt(dot_prod_self(m_w));
    double residual = 0.0, residual_prev;

    if(utils::fpu::is_fpu_error(rp_norm))
    {
        fprintf(stderr, "Error: FPU error detected in right part!\n");
        m_precond->solve_Q(m_x0, solution);
        delete [] m_x0;
        return;
    }

    bool finished = false;

    std::size_t iter;
    for(iter = 0; iter <= max_iter && !finished; iter++)
    {
        discr = sqrt(dot_prod_self(m_r));
        if(utils::fpu::is_fpu_error(discr))
        {
            fprintf(stderr, "Error: FPU error detected in (r, r)!\n");
            m_precond->solve_Q(m_x0, solution);
            delete [] m_x0;
            return;
        }

        residual_prev = residual;
        residual = discr / rp_norm;
        //if(iter%10 == 0)
        {
            printf("GMRES_Complex<%s>(%lu) Residual:\t%5lu\t%.3e\r", precond_name, static_cast<unsigned long>(m_m), static_cast<unsigned long>(iter), residual);
            fflush(stdout);
        }

        if(residual > eps && fabs(residual - residual_prev) / (residual) > 1e-15)
        {
            // V
            for(std::size_t i = 0; i < m_n; i++)
                m_VT[0][i] = m_r[i] / discr;

            // d
            m_d[0] = discr;
            for(std::size_t i = 1; i <= m_m_curr; i++)
                m_d[i] = 0.0;

            for(std::size_t mu = 0; mu < m_m_curr; mu++)
            {
                // w
                solve_QAS(m_VT[mu], m_w, m_t);

                // H_lambda_mu
                for(std::size_t lambda = 0; lambda <= mu; lambda++)
                    m_H[lambda][mu] = dot_prod(m_VT[lambda], m_w);

                // tilde_v_mu+1
                for(std::size_t lambda = 0; lambda <= mu; lambda++)
                    for(std::size_t i = 0; i < m_n; i++)
                        m_w[i] -= m_VT[lambda][i] * m_H[lambda][mu];

                // H_mu+1_mu
                double vnorm = sqrt(dot_prod_self(m_w));
                m_H[mu + 1][mu] = vnorm;

                if(std::abs(vnorm) < 1e-15)
                {
                    if(mu > 1)  m_m_curr = mu;
                    else        m_m_curr = 1;
                    break;
                }

                // v_mu+1
                if(mu + 1 < m_m_curr)
                    for(std::size_t i = 0; i < m_n; i++)
                        m_VT[mu + 1][i] = m_w[i] / vnorm;
            }

            // H(i), d(i)
            for(std::size_t i = 0; i < m_m_curr; i++)
            {
                double t1 = m_H[i + 1][i].real();
                double t2 = abs(m_H[i][i]);
                double t3 = 1.0 / sqrt(t2 * t2 + t1 * t1);
                double s = t1 * t3;
                std::complex<double> c = m_H[i][i] * t3;

                for(std::size_t j = i; j < m_m_curr; j++)
                {
                    std::complex<double> H_ij = m_H[i][j];
                    m_H[i][j] = H_ij * conj(c) + m_H[i + 1][j] * s;
                    m_H[i + 1][j] = - H_ij * s + m_H[i + 1][j] * c;
                }
                std::complex<double> d_i = m_d[i];
                m_d[i] = d_i * conj(c) + m_d[i + 1] * s;
                m_d[i + 1] = - d_i * s + m_d[i + 1] * c;
            }

            // Hz = d
            for(std::size_t i = 0; i < m_m_curr; i++)
            {
                m_d[i] = m_d[i] / m_H[i][i];
                for(std::size_t j = i + 1; j < m_m_curr; j++)
                    m_H[i][j] /= m_H[i][i];
            }
            for(std::size_t i = m_m_curr - 1; i > 0; i--)
                for(std::size_t j = i - 1, ju = i; ju > 0; j--, ju--)
                    m_d[j] -= m_H[j][i] * m_d[i];

            // x = x + Vz
            for(std::size_t i = 0; i < m_n; i++)
                for(std::size_t j = 0; j < m_m_curr; j++)
                    m_x0[i] += m_VT[j][i] * m_d[j];

            // r
            copy_vec(m_x0, m_t);
            m_precond->solve_Q(m_t, m_t);
            mul_matrix(m_t, m_r);
            for(std::size_t i = 0; i < m_n; i++)
                m_t[i] = m_rp[i] - m_r[i];
            m_precond->solve_S(m_t, m_r);
        }
        else
            finished = true;
    }

    m_precond->solve_Q(m_x0, solution);
    delete [] m_x0;

//    mul_matrix(solution, m_r);
//    for(std::size_t i = 0; i < m_n; i++)
//        m_r[i] = m_rp[i] - m_r[i];
//    discr = sqrt(dot_prod_self(m_r));
//    rp_norm = sqrt(dot_prod_self(m_rp));
    printf("GMRES_Complex<%s>(%lu) Residual:\t%5lu\t%.3e\n", precond_name_str.c_str(), static_cast<unsigned long>(m_m), static_cast<unsigned long>(iter) - 1, discr / rp_norm);

    if(iter >= max_iter)
        printf("Soulution can`t found, iteration limit exceeded!\n");
    else if(residual > eps)
        printf("Soulution can`t found, stagnation detected!\n");
}

GMRES_Complex::GMRES_Complex()
{
    m_d = m_r = m_w = m_t = m_x0 = NULL;
    m_H = m_VT = NULL;
    m_m = m_m_curr = 0;
}

GMRES_Complex::~GMRES_Complex()
{
    delete [] m_r;
    delete [] m_w;
    delete [] m_t;
    delete [] m_d;
    if(m_H && m_m)
        for(std::size_t i = 0; i <= m_m; i++)
            delete [] m_H[i];
    delete [] m_H;
    if(m_VT && m_m)
        for(std::size_t i = 0; i < m_m; i++)
            delete [] m_VT[i];
    delete [] m_VT;
}

}}}}} // namespace fem_core::solvers::CSLR::symmetric::complex
