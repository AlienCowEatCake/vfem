#if defined(_MSC_VER) && !defined(_CRT_SECURE_NO_WARNINGS)
#define _CRT_SECURE_NO_WARNINGS
#endif
#include "COCR_Smooth.h"
#include <cstdio>
#include <cmath>
#include "../../../../../utils/fpu.h"

namespace fem_core { namespace solvers { namespace CSLR { namespace symmetric { namespace complex {

void COCR_Smooth::init(const std::size_t * gi, const std::size_t * gj, const std::complex<double> * di,
                       const std::complex<double> * gg, std::size_t n,
                       preconditioners::preconditioner_interface<std::complex<double>, std::size_t> * precond)
{
    COCR::init(gi, gj, di, gg, n, precond);
    delete [] m_xs;
    delete [] m_rs;
    m_xs = new std::complex<double> [m_n];
    m_rs = new std::complex<double> [m_n];
}

double COCR_Smooth::dot_prod_real(const std::complex<double> * a, const std::complex<double> * b) const
{
    double d_p = 0.0;
    for(std::size_t i = 0; i < m_n; i++)
        d_p += a[i].real() * b[i].real() + a[i].imag() * b[i].imag();
    return d_p;
}

void COCR_Smooth::solve(std::complex<double> * solution, const std::complex<double> * rp_s,
                        double eps, std::size_t max_iter)
{
    m_rp = rp_s;
    std::string precond_name_str = m_precond->get_name();
    const char * precond_name = precond_name_str.c_str();

    m_x0 = solution;
    for(std::size_t i = 0; i < m_n; i++)
        m_xs[i] = m_x0[i];

    mul_matrix(m_x0, m_r);

    for(std::size_t i = 0; i < m_n ; i++)
        m_rs[i] = m_r[i] = m_rp[i] - m_r[i];

    solve_SQ(m_r, m_s);
    mul_matrix(m_s, m_z);
    solve_SQ(m_z, m_w);
    for(std::size_t i = 0; i < m_n; i++)
    {
        m_p[i] = m_s[i];
        m_a[i] = m_z[i];
    }
    std::complex<double> dp_as = dot_prod_nocj(m_a, m_s), dp_as_new;
    std::complex<double> alpha, beta;
    double discr = 0.0, rp_norm, residual = 0.0, /*residual_old,*/ eta;

    rp_norm = sqrt(dot_prod_self(m_rp));
    if(utils::fpu::is_fpu_error(rp_norm))
    {
        fprintf(stderr, "Error: FPU error detected in right part!\n");
        return;
    }

    bool finished = false;
    std::size_t iter;
    for(iter = 0; iter <= max_iter && !finished; iter++)
    {
        discr = sqrt(dot_prod_self(m_rs));
        if(utils::fpu::is_fpu_error(discr))
        {
            fprintf(stderr, "Error: FPU error detected in (r, r)!\n");
            for(std::size_t i = 0; i < m_n; i++)
                solution[i] = m_xs[i];
            return;
        }

        /*residual_old = residual;*/
        residual = discr / rp_norm;
        //if(iter%10 == 0)
        {
            printf("COCR_Smooth<%s> Residual:\t%5lu\t%.3e\r", precond_name, static_cast<unsigned long>(iter), residual);
            fflush(stdout);
        }

        if(residual > eps/* && fabs(residual - residual_old) / (residual) > 1e-15*/)
        {
            alpha = dp_as / dot_prod_nocj(m_z, m_w);
            for(std::size_t i = 0; i < m_n ; i++)
            {
                m_x0[i] += alpha * m_p[i];
                m_r[i] -= alpha * m_z[i];
                m_s[i] -= alpha * m_w[i];
                m_w[i] = m_r[i] - m_rs[i]; // r - s, сглаживатель
            }
            eta = - dot_prod_real(m_rs, m_w) / dot_prod_self(m_w);
            if(eta < 0.0) eta = 0.0;
            else if(eta > 1.0) eta = 1.0;
            for(std::size_t i = 0; i < m_n ; i++)
            {
                double eta1 = 1.0 - eta;
                m_xs[i] = eta1 * m_xs[i] + eta * m_x0[i];
                m_rs[i] = eta1 * m_rs[i] + eta * m_r[i];
            }
            mul_matrix(m_s, m_a);
            dp_as_new = dot_prod_nocj(m_a, m_s);
            beta = dp_as_new / dp_as;
            dp_as = dp_as_new;
            for(std::size_t i = 0; i < m_n ; i++)
            {
                m_p[i] = m_s[i] + beta * m_p[i];
                m_z[i] = m_a[i] + beta * m_z[i];
            }
            solve_SQ(m_z, m_w);
        }
        else
            finished = true;
    }

//    mul_matrix(m_xs, m_r);
//    for(std::size_t i = 0; i < m_n; i++)
//        m_r[i] = m_rp[i] - m_r[i];
//    discr = sqrt(dot_prod_self(m_r));
    printf("COCR_Smooth<%s> Residual:\t%5lu\t%.3e\n", precond_name_str.c_str(), static_cast<unsigned long>(iter) - 1, discr / rp_norm);

    if(iter >= max_iter)
        printf("Soulution can`t found, iteration limit exceeded!\n");
    else if(residual > eps)
        printf("Soulution can`t found, stagnation detected!\n");

    for(std::size_t i = 0; i < m_n; i++)
        solution[i] = m_xs[i];
}

COCR_Smooth::COCR_Smooth()
{
    m_xs = m_rs = NULL;
}

COCR_Smooth::~COCR_Smooth()
{
    delete [] m_xs;
    delete [] m_rs;
}

}}}}} // namespace fem_core::solvers::CSLR::symmetric::complex
