#if defined(_MSC_VER) && !defined(_CRT_SECURE_NO_WARNINGS)
#define _CRT_SECURE_NO_WARNINGS
#endif
#include "BiCGStab_Complex.h"
#include <cstdio>
#include <cmath>
#include "../../../../../utils/fpu.h"

namespace core { namespace solvers { namespace CSLR { namespace symmetric { namespace complex {

BiCGStab_Complex::BiCGStab_Complex()
{
    m_r = m_v = m_s = m_p = m_t = m_r2 = NULL;
}

BiCGStab_Complex::~BiCGStab_Complex()
{
    delete [] m_r;
    delete [] m_v;
    delete [] m_s;
    delete [] m_p;
    delete [] m_t;
    delete [] m_r2;
}

void BiCGStab_Complex::init(const std::size_t * gi, const std::size_t * gj, const std::complex<double> * di,
                            const std::complex<double> * gg, std::size_t n,
                            preconditioners::preconditioner_interface<std::complex<double>, std::size_t> * precond)
{
    /// @todo Добавить работу с предобуславливателем
    (void)(precond);
    m_gi = gi;
    m_gj = gj;
    m_di = di;
    m_gg = gg;
    m_n = n;

    delete [] m_r;
    delete [] m_v;
    delete [] m_s;
    delete [] m_p;
    delete [] m_t;
    delete [] m_r2;

    m_r = new std::complex<double> [m_n];
    m_v = new std::complex<double> [m_n];
    m_s = new std::complex<double> [m_n];
    m_p = new std::complex<double> [m_n];
    m_t = new std::complex<double> [m_n];
    m_r2 = new std::complex<double> [m_n];
}

std::complex<double> BiCGStab_Complex::dot_prod(const std::complex<double> * a, const std::complex<double> * b) const
{
    std::complex<double> res;
    for(std::size_t i = 0; i < m_n; i++)
    {
        res += conj(a[i]) * b[i];
    }
    return res;
}

std::complex<double> BiCGStab_Complex::dot_prod_nocj(const std::complex<double> * a, const std::complex<double> * b) const
{
    std::complex<double> res;
    for(std::size_t i = 0; i < m_n; i++)
    {
        res += a[i] * b[i];
    }
    return res;
}

double BiCGStab_Complex::dot_prod_self(const std::complex<double> * a) const
{
    double res = 0.0;
    for(std::size_t i = 0; i < m_n; i++)
    {
        double re = a[i].real();
        double im = a[i].imag();
        res += re * re + im * im;
    }
    return res;
}

void BiCGStab_Complex::mul_matrix(const std::complex<double> * f, std::complex<double> * x) const
{
    for(std::size_t i = 0; i < m_n; i++)
    {
        std::complex<double> v_el = f[i];
        x[i] = m_di[i] * v_el;
        for(std::size_t k = m_gi[i], k1 = m_gi[i + 1]; k < k1; k++)
        {
            std::size_t j = m_gj[k];
            x[i] += m_gg[k] * f[j];
            x[j] += m_gg[k] * v_el;
        }
    }
}

void BiCGStab_Complex::solve(std::complex<double> * solution, const std::complex<double> * rp,
                             double gamma, std::size_t max_iter)
{
    double eps = gamma * gamma;

    std::complex<double> omega, alpha, ro, beta;
    double rp_norm = dot_prod_self(rp);
    if(utils::fpu::is_fpu_error(rp_norm))
    {
        fprintf(stderr, "Error: FPU error detected in right part!\n");
        return;
    }
    m_x0 = solution;

    mul_matrix(solution, m_t);
    for(std::size_t i = 0; i < m_n; i++)
    {
        m_p[i] = m_r[i] = rp[i] - m_t[i];
        m_r2[i] = 1.0;
        m_v[i] = 0;
    }
    ro = alpha = omega = 1.0;

    std::complex<double> a1 = 0.0, a2 = 0.0, a3;

    bool not_end = true;
    double discr = 0.0;
    std::size_t iter;
    for(iter = 0; iter < max_iter && not_end; iter++)
    {
        discr = dot_prod_self(m_r);
        if(utils::fpu::is_fpu_error(discr))
        {
            fprintf(stderr, "Error: FPU error detected in (r, r)!\n");
            return;
        }

        //if(iter%10 == 0)
        {
            printf("BiCGStab_Complex Residual:\t%5lu\t%.3e\r", static_cast<unsigned long>(iter), sqrt(discr / rp_norm));
            fflush(stdout);
        }

        if(discr / rp_norm > eps)
        {
            mul_matrix(m_p, m_v);
            a1 = dot_prod_nocj(m_r, m_r2);
            a2 = dot_prod_nocj(m_v, m_r2);
            alpha = a1 / a2;
            for(std::size_t i = 0; i < m_n; i++)
                m_s[i] = m_r[i] - alpha * m_v[i];

            mul_matrix(m_s, m_t);
            a2 = dot_prod_nocj(m_s, m_t);
            a3 = dot_prod_nocj(m_t, m_t);
            omega = a2 / a3;

            for(std::size_t i = 0; i < m_n; i++)
            {
                m_x0[i] += alpha * m_p[i] + omega * m_s[i];
                m_r[i] = m_s[i] - omega * m_t[i];
            }

            a2 = dot_prod_nocj(m_r, m_r2);
            beta = a2 / a1 * alpha / omega;

            for(std::size_t i = 0; i < m_n; i++)
                m_p[i] = m_r[i] + beta * (m_p[i] - omega * m_v[i]);
        }
        else
            not_end = false;
    }

//    mul_matrix(m_x0, m_r);
//    for(std::size_t i = 0; i < m_n; i++)
//        m_r[i] = rp[i] - m_r[i];
//    discr = dot_prod_self(m_r);
    printf("BiCGStab_Complex Residual:\t%5lu\t%.3e\n", static_cast<unsigned long>(iter) - 1, sqrt(discr / rp_norm));

    if(iter >= max_iter)
        printf("Soulution can`t found, iteration limit exceeded!\n");
}

}}}}} // namespace core::solvers::CSLR::symmetric::complex
