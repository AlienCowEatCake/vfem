#if defined(_MSC_VER) && !defined(_CRT_SECURE_NO_WARNINGS)
#define _CRT_SECURE_NO_WARNINGS
#endif
#include "BiCG_Complex_Smooth.h"
#include <cstdio>
#include <cmath>
#include "../../../../../utils/fpu.h"

namespace core { namespace solvers { namespace CSRC { namespace symmetric { namespace complex {

BiCG_Complex_Smooth::BiCG_Complex_Smooth()
{
    m_xs = m_rs = NULL;
}

BiCG_Complex_Smooth::~BiCG_Complex_Smooth()
{
    delete [] m_xs;
    delete [] m_rs;
}

void BiCG_Complex_Smooth::init(const std::size_t * gi, const std::size_t * gj, const std::complex<double> * di,
                               const std::complex<double> * gg, std::size_t n,
                               preconditioners::preconditioner_interface<std::complex<double>, std::size_t> * precond)
{
    /// @todo Добавить работу с предобуславливателем
    (void)(precond);
    BiCG_Complex::init(gi, gj, di, gg, n, precond);
    delete [] m_xs;
    delete [] m_rs;
    m_xs = new std::complex<double> [m_n];
    m_rs = new std::complex<double> [m_n];
}

double BiCG_Complex_Smooth::dot_prod_real(const std::complex<double> * a, const std::complex<double> * b) const
{
    double d_p = 0.0;
    for(std::size_t i = 0; i < m_n; i++)
        d_p += a[i].real() * b[i].real() + a[i].imag() * b[i].imag();
    return d_p;
}

void BiCG_Complex_Smooth::solve(std::complex<double> * solution, const std::complex<double> * rp, double gamma, std::size_t max_iter)
{
    double eps = gamma * gamma;

    std::complex<double> dp1, dp2, alpha, beta;
    double rp_norm = dot_prod_self(rp);
    m_x0 = solution;

    mul_matrix(solution, m_t);
    for(std::size_t i = 0; i < m_n; i++)
    {
        m_xs[i] = solution[i];
        m_rs[i] = m_r[i] = rp[i] - m_t[i];
        m_p[i] = m_z[i] = m_s[i] = m_r[i];
    }

    bool not_end = true;
    double discr = 0.0;
    std::size_t iter;

    dp1 = dot_prod_nocj(m_p, m_r);
    if(utils::fpu::is_fpu_error(dp1.real()))
    {
        fprintf(stderr, "Error: FPU error detected in right part!\n");
        return;
    }

    for(iter = 0; iter < max_iter && not_end; iter++)
    {
        discr = dot_prod_self(m_rs);
        if(utils::fpu::is_fpu_error(discr))
        {
            fprintf(stderr, "Error: FPU error detected in (r, r)!\n");
            for(std::size_t i = 0; i < m_n; i++)
                solution[i] = m_xs[i];
            return;
        }

        //if(iter%10 == 0)
        {
            printf("BiCG_Complex_Smooth Residual:\t%5lu\t%.3e\r", (unsigned long)iter, sqrt(discr / rp_norm));
            fflush(stdout);
        }

        if(discr / rp_norm > eps)
        {
            mul_matrix(m_z, m_t);

            alpha = dp1 / dot_prod_nocj(m_s, m_t);

            mul_matrix(m_s, m_t1);

            for(std::size_t i = 0; i < m_n; i++)
            {
                m_x0[i] += alpha * m_z[i];
                m_r[i] -= alpha * m_t[i];
                m_p[i] -= alpha * m_t1[i];
                m_t[i] = m_r[i] - m_rs[i]; // r - s, сглаживатель
            }

            double eta = - dot_prod_real(m_rs, m_t) / dot_prod_self(m_t);
            if(eta < 0.0) eta = 0.0;
            else if(eta > 1.0) eta = 1.0;
            for(std::size_t i = 0; i < m_n ; i++)
            {
                double eta1 = 1.0 - eta;
                m_xs[i] = eta1 * m_xs[i] + eta * m_x0[i];
                m_rs[i] = eta1 * m_rs[i] + eta * m_r[i];
            }

            dp2 = dot_prod_nocj(m_p, m_r);
            beta = dp2 / dp1;
            dp1 = dp2;

            for(std::size_t i = 0; i < m_n; i++)
            {
                m_z[i] = m_r[i] + beta * m_z[i];
                m_s[i] = m_p[i] + beta * m_s[i];
            }
        }
        else
            not_end = false;
    }

//    mul_matrix(m_xs, m_r);
//    for(std::size_t i = 0; i < m_n; i++)
//        m_r[i] = rp[i] - m_r[i];
//    discr = dot_prod_self(m_r);
    printf("BiCG_Complex_Smooth Residual:\t%5lu\t%.3e\n", (unsigned long)iter - 1, sqrt(discr / rp_norm));

    if(iter >= max_iter)
        printf("Soulution can`t found, iteration limit exceeded!\n");

    for(std::size_t i = 0; i < m_n; i++)
        solution[i] = m_xs[i];
}

}}}}} // namespace core::solvers::CSRC::symmetric::complex
