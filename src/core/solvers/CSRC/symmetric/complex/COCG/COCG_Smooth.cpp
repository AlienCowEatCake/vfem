#if defined(_MSC_VER) && !defined(_CRT_SECURE_NO_WARNINGS)
#define _CRT_SECURE_NO_WARNINGS
#endif
#include "COCG_Smooth.h"
#include <cstdio>
#include <cmath>
#include "../../../../../utils/fpu.h"

namespace core { namespace solvers { namespace CSRC { namespace symmetric { namespace complex {

void COCG_Smooth::init(const std::size_t * gi, const std::size_t * gj, const std::complex<double> * di,
                       const std::complex<double> * gg, std::size_t n,
                       preconditioners::preconditioner_interface<std::complex<double>, std::size_t> * precond)
{
    COCG::init(gi, gj, di, gg, n, precond);
    delete [] m_xs;
    delete [] m_rs;
    m_xs = new std::complex<double> [m_n];
    m_rs = new std::complex<double> [m_n];
}

double COCG_Smooth::dot_prod_real(const std::complex<double> * a, const std::complex<double> * b) const
{
    double d_p = 0.0;
    for(std::size_t i = 0; i < m_n; i++)
        d_p += a[i].real() * b[i].real() + a[i].imag() * b[i].imag();
    return d_p;
}

void COCG_Smooth::solve(std::complex<double> * solution, const std::complex<double> * rp,
                        double eps, std::size_t max_iter)
{
    eps *= eps;
    std::string precond_name_str = m_precond->get_name();
    const char * precond_name = precond_name_str.c_str();

    m_rp = rp;

    m_x0 = solution;
    for(std::size_t i = 0; i < m_n; i++)
        m_xs[i] = m_x0[i];

    mul_matrix(m_x0, m_r);

    for(std::size_t i = 0; i < m_n ; i++)
        m_rs[i] = m_r[i] = m_rp[i] - m_r[i];

    solve_SQ(m_r, m_z);
    for(std::size_t i = 0; i < m_n; i++)
        m_p[i] = m_z[i];

    std::complex<double> alpha, beta, prod_1, prod_2;
    double discr = 0.0, rp_norm, eta;

    rp_norm = dot_prod_self(m_rp);
    if(utils::fpu::is_fpu_error(rp_norm))
    {
        fprintf(stderr, "Error: FPU error detected in right part!\n");
        return;
    }

    prod_1 = dot_prod_nocj(m_p, m_r);

    bool finished = false;

    std::size_t iter;
    for(iter = 0; iter <= max_iter && !finished; iter++)
    {
        discr = dot_prod_self(m_rs);
        if(utils::fpu::is_fpu_error(discr))
        {
            fprintf(stderr, "Error: FPU error detected in (r, r)!\n");
            for(std::size_t i = 0; i < m_n; i++)
                solution[i] = m_xs[i];
            return;
        }

        double residual = discr / rp_norm;
//        if(iter%10 == 0)
        {
            printf("COCG_Smooth<%s> Residual:\t%5lu\t%.3e\r", precond_name, (unsigned long)iter, sqrt(residual));
            fflush(stdout);
        }

        if(residual > eps)
        {
            mul_matrix(m_z, m_s);

            alpha = prod_1 / dot_prod_nocj(m_s, m_z);

            for(std::size_t i = 0; i < m_n ; i++)
            {
                m_x0[i] += alpha * m_z[i];
                m_r[i] -= alpha * m_s[i];
                m_s[i] = m_r[i] - m_rs[i]; // r - s, сглаживатель
            }

            eta = - dot_prod_real(m_rs, m_s) / dot_prod_self(m_s);
            if(eta < 0.0) eta = 0.0;
            else if(eta > 1.0) eta = 1.0;
            for(std::size_t i = 0; i < m_n ; i++)
            {
                double eta1 = 1.0 - eta;
                m_xs[i] = eta1 * m_xs[i] + eta * m_x0[i];
                m_rs[i] = eta1 * m_rs[i] + eta * m_r[i];
            }

            solve_SQ(m_r, m_p);
            prod_2 = dot_prod_nocj(m_p, m_r);

            beta = prod_2 / prod_1;

            prod_1 = prod_2;

            for(std::size_t i = 0; i < m_n; i++)
                m_z[i] = m_p[i] + beta * m_z[i];
        }
        else
            finished = true;
    }

//    mul_matrix(m_xs, m_r);
//    for(std::size_t i = 0; i < m_n; i++)
//        m_r[i] = m_rp[i] - m_r[i];
//    discr = dot_prod_self(m_r);
    printf("COCG_Smooth<%s> Residual:\t%5lu\t%.3e\n", precond_name_str.c_str(), (unsigned long)iter - 1, sqrt(discr / rp_norm));

    if(iter >= max_iter)
        printf("Soulution can`t found, iteration limit exceeded!\n");

    for(std::size_t i = 0; i < m_n; i++)
        solution[i] = m_xs[i];
}

COCG_Smooth::COCG_Smooth()
{
    m_xs = m_rs = NULL;
}

COCG_Smooth::~COCG_Smooth()
{
    delete [] m_xs;
    delete [] m_rs;
}

}}}}} // namespace core::solvers::CSRC::symmetric::complex
