#if defined(_MSC_VER) && !defined(_CRT_SECURE_NO_WARNINGS)
#define _CRT_SECURE_NO_WARNINGS
#endif
#include "COCG.h"
#include <cstdio>
#include <cmath>
#include "../../../../../utils/fpu.h"

namespace core { namespace solvers { namespace CSRC { namespace symmetric { namespace complex {

void COCG::init(const std::size_t * gi, const std::size_t * gj, const std::complex<double> * di,
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
    delete [] m_z;
    delete [] m_p;
    delete [] m_s;

    m_r = new std::complex<double> [m_n];
    m_z = new std::complex<double> [m_n];
    m_p = new std::complex<double> [m_n];
    m_s = new std::complex<double> [m_n];
}

std::complex<double> COCG::dot_prod_nocj(const std::complex<double> * a, const std::complex<double> * b) const
{
    std::complex<double> d_p = 0.0;
    for(std::size_t i = 0; i < m_n; i++)
        d_p += a[i] * b[i];
    return d_p;
}

double COCG::dot_prod_self(const std::complex<double> * a) const
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

void COCG::mul_matrix(const std::complex<double> * f, std::complex<double> * x) const
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

void COCG::solve_SQ(const std::complex<double> * f, std::complex<double> * x) const
{
    m_precond->solve_S(f, x);
    m_precond->solve_Q(x, x);
}

void COCG::solve(std::complex<double> * solution, const std::complex<double> * rp,
                 double eps, std::size_t max_iter)
{
    eps *= eps;
    std::string precond_name_str = m_precond->get_name();
    const char * precond_name = precond_name_str.c_str();

    m_rp = rp;

    m_x0 = solution;

    mul_matrix(m_x0, m_r);

    for(std::size_t i = 0; i < m_n ; i++)
        m_r[i] = m_rp[i] - m_r[i];

    solve_SQ(m_r, m_z);
    for(std::size_t i = 0; i < m_n; i++)
        m_p[i] = m_z[i];

    std::complex<double> alpha, beta, prod_1, prod_2;
    double discr = 0.0, rp_norm;

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
        discr = dot_prod_self(m_r);
        if(utils::fpu::is_fpu_error(discr))
        {
            fprintf(stderr, "Error: FPU error detected in (r, r)!\n");
            return;
        }

        double residual = discr / rp_norm;
//        if(iter%10 == 0)
        {
            printf("COCG<%s> Residual:\t%5lu\t%.3e\r", precond_name, (unsigned long)iter, sqrt(residual));
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

//    mul_matrix(m_x0, m_r);
//    for(std::size_t i = 0; i < m_n; i++)
//        m_r[i] = m_rp[i] - m_r[i];
//    discr = dot_prod_self(m_r);
    printf("COCG<%s> Residual:\t%5lu\t%.3e\n", precond_name_str.c_str(), (unsigned long)iter - 1, sqrt(discr / rp_norm));

    if(iter >= max_iter)
        printf("Soulution can`t found, iteration limit exceeded!\n");
}

COCG::COCG()
{
    m_r = m_x0 = m_z = m_p = m_s = NULL;
}

COCG::~COCG()
{
    delete [] m_r;
    delete [] m_z;
    delete [] m_p;
    delete [] m_s;
}

}}}}} // namespace core::solvers::CSRC::symmetric::complex
