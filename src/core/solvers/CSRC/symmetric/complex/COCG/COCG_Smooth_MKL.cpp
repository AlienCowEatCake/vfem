#if defined(_MSC_VER) && !defined(_CRT_SECURE_NO_WARNINGS)
#define _CRT_SECURE_NO_WARNINGS
#endif
#include "COCG_Smooth_MKL.h"
#include <cstdio>
#include <cmath>
#include <algorithm>
#include "../../../../../utils/fpu.h"
#include "../../../../../wrappers/omp_wrapper.h"
#include "../../../../../wrappers/mkl_wrapper.h"

namespace core { namespace solvers { namespace CSRC { namespace symmetric { namespace complex {

typedef wrappers::omp::omp_int omp_int;

// y := x
inline void normal_zcopy(std::size_t n, const std::complex<double> * x, std::complex<double> * y)
{
    MKL_INT m = static_cast<MKL_INT>(n);
    MKL_Complex16 * x1 = reinterpret_cast<MKL_Complex16 *>(const_cast<std::complex<double> *>(x));
    MKL_Complex16 * y1 = reinterpret_cast<MKL_Complex16 *>(const_cast<std::complex<double> *>(y));
    cblas_zcopy(m, x1, 1, y1, 1);
}

// y := a * x + y
inline void normal_zaxpy(std::size_t n, const std::complex<double> & a, const std::complex<double> * x, std::complex<double> * y)
{
    MKL_INT m = static_cast<MKL_INT>(n);
    MKL_Complex16 * a1 = reinterpret_cast<MKL_Complex16 *>(const_cast<std::complex<double> *>(&a));
    MKL_Complex16 * x1 = reinterpret_cast<MKL_Complex16 *>(const_cast<std::complex<double> *>(x));
    MKL_Complex16 * y1 = reinterpret_cast<MKL_Complex16 *>(const_cast<std::complex<double> *>(y));
    cblas_zaxpy(m, a1, x1, 1, y1, 1);
}

void COCG_Smooth_MKL::init(const std::size_t * gi, const std::size_t * gj, const std::complex<double> * di,
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
    delete [] m_xs;
    delete [] m_rs;

    m_r = new std::complex<double> [m_n];
    m_z = new std::complex<double> [m_n];
    m_p = new std::complex<double> [m_n];
    m_s = new std::complex<double> [m_n];
    m_xs = new std::complex<double> [m_n];
    m_rs = new std::complex<double> [m_n];

    m_m = static_cast<MKL_INT>(m_n);
    m_ia = new MKL_INT[m_n + 1];
    m_ia[0] = 1;
    for(MKL_INT i = 1; i <= m_m; i++)
        m_ia[i] = m_ia[i - 1] + static_cast<MKL_INT>(m_gi[i] - m_gi[i - 1]) + 1;
    m_ja = new MKL_INT[m_ia[m_n]];
    m_aa = new MKL_Complex16[m_ia[m_n]];
#pragma omp parallel for
    for(omp_int i = 0; i < static_cast<omp_int>(m_n); i++)
    {
        std::size_t i_b = m_gi[i], i_e = m_gi[i + 1];
        for(std::size_t j = i_b; j < i_e; j++)
        {
            std::size_t ind = static_cast<std::size_t>(m_ia[i]) + (j - m_gi[i]) - 1;
            m_ja[ind] = static_cast<MKL_INT>(m_gj[j]) + 1;
            m_aa[ind].real = m_gg[j].real();
            m_aa[ind].imag = m_gg[j].imag();
        }
        std::size_t ind = static_cast<std::size_t>(m_ia[i + 1] - 2);
        m_ja[ind] = static_cast<MKL_INT>(i) + 1;
        m_aa[ind].real = m_di[i].real();
        m_aa[ind].imag = m_di[i].imag();
    }
}

std::complex<double> COCG_Smooth_MKL::dot_prod_nocj(const std::complex<double> * a, const std::complex<double> * b) const
{
    std::complex<double> result;
    MKL_Complex16 * a1 = reinterpret_cast<MKL_Complex16 *>(const_cast<std::complex<double> *>(a));
    MKL_Complex16 * b1 = reinterpret_cast<MKL_Complex16 *>(const_cast<std::complex<double> *>(b));
    MKL_Complex16 * pres = reinterpret_cast<MKL_Complex16 *>(&result);
    cblas_zdotu_sub(m_m, a1, 1, b1, 1, pres);
    return result;
}

double COCG_Smooth_MKL::dot_prod_self(const std::complex<double> * a) const
{
    MKL_Complex16 * a1 = reinterpret_cast<MKL_Complex16 *>(const_cast<std::complex<double> *>(a));
    double result = cblas_dznrm2(m_m, a1, 1);
    return result * result;
}

double COCG_Smooth_MKL::dot_prod_real(const std::complex<double> * a, const std::complex<double> * b) const
{
    double d_p = 0.0;
#pragma omp parallel for reduction(+ : d_p)
    for(omp_int i = 0; i < static_cast<omp_int>(m_n); i++)
        d_p += a[i].real() * b[i].real() + a[i].imag() * b[i].imag();
    return d_p;
}

void COCG_Smooth_MKL::mul_matrix(const std::complex<double> * f, std::complex<double> * x) const
{
    MKL_Complex16 * in_v = reinterpret_cast<MKL_Complex16 *>(const_cast<std::complex<double> *>(f));
    MKL_Complex16 * out_v = reinterpret_cast<MKL_Complex16 *>(const_cast<std::complex<double> *>(x));
    mkl_zcsrsymv("L", &m_m, m_aa, m_ia, m_ja, in_v, out_v);
}

void COCG_Smooth_MKL::solve_SQ(const std::complex<double> * f, std::complex<double> * x) const
{
    m_precond->solve_S(f, x);
    m_precond->solve_Q(x, x);
}

void COCG_Smooth_MKL::solve(std::complex<double> * solution, const std::complex<double> * rp,
                               double eps, std::size_t max_iter)
{
    eps *= eps;
    std::string precond_name_str = m_precond->get_name();
    const char * precond_name = precond_name_str.c_str();

    m_rp = rp;

    m_x0 = solution;
    mul_matrix(m_x0, m_r);

#pragma omp parallel for
    for(omp_int i = 0; i < static_cast<omp_int>(m_n); i++)
    {
        m_r[i] = m_rp[i] - m_r[i];
    }
    solve_SQ(m_r, m_z);
    normal_zcopy(m_n, solution, m_xs);
    normal_zcopy(m_n, m_r, m_rs);
    normal_zcopy(m_n, m_z, m_p);

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
            normal_zcopy(m_n, m_xs, solution);
            return;
        }

        double residual = discr / rp_norm;
//        if(iter%10 == 0)
        {
            printf("COCG_Smooth_MKL<%s> [%d] Residual:\t%5lu\t%.3e\r", precond_name, m_num_threads, (unsigned long)iter, sqrt(residual));
            fflush(stdout);
        }

        if(residual > eps)
        {
            mul_matrix(m_z, m_s);

            alpha = prod_1 / dot_prod_nocj(m_s, m_z);

            normal_zaxpy(m_n, alpha, m_z, m_x0);
            normal_zaxpy(m_n, -alpha, m_s, m_r);
#pragma omp parallel for
            for(omp_int i = 0; i < static_cast<omp_int>(m_n); i++)
            {
                m_s[i] = m_r[i] - m_rs[i]; // r - s, сглаживатель
            }

            eta = - dot_prod_real(m_rs, m_s) / dot_prod_self(m_s);
            if(eta < 0.0) eta = 0.0;
            else if(eta > 1.0) eta = 1.0;
#pragma omp parallel for
            for(omp_int i = 0; i < static_cast<omp_int>(m_n); i++)
            {
                double eta1 = 1.0 - eta;
                m_xs[i] = eta1 * m_xs[i] + eta * m_x0[i];
                m_rs[i] = eta1 * m_rs[i] + eta * m_r[i];
            }

            solve_SQ(m_r, m_p);
            prod_2 = dot_prod_nocj(m_p, m_r);

            beta = prod_2 / prod_1;

            prod_1 = prod_2;

            normal_zaxpy(m_n, beta, m_z, m_p);
            std::swap(m_z, m_p);
        }
        else
            finished = true;
    }

//    mul_matrix(m_xs, m_r);
//#pragma omp parallel for
//    for(omp_int i = 0; i < static_cast<omp_int>(m_n); i++)
//        m_r[i] = m_rp[i] - m_r[i];
//    discr = dot_prod_self(m_r);
    printf("COCG_Smooth_MKL<%s> [%d] Residual:\t%5lu\t%.3e\n", precond_name_str.c_str(), m_num_threads, (unsigned long)iter - 1, sqrt(discr / rp_norm));

    if(iter >= max_iter)
        printf("Soulution can`t found, iteration limit exceeded!\n");

    normal_zcopy(m_n, m_xs, solution);
}

COCG_Smooth_MKL::COCG_Smooth_MKL()
{
    wrappers::mkl::wrapper_mkl_set_env_max_threads();
    m_num_threads = std::max(mkl_get_max_threads(), omp_get_max_threads());
    m_r = m_x0 = m_z = m_p = m_s = m_xs = m_rs = NULL;
    m_ia = m_ja = NULL;
    m_aa = NULL;
}

COCG_Smooth_MKL::~COCG_Smooth_MKL()
{
    delete [] m_r;
    delete [] m_z;
    delete [] m_p;
    delete [] m_s;
    delete [] m_xs;
    delete [] m_rs;
    delete [] m_ia;
    delete [] m_ja;
    delete [] m_aa;
}

}}}}} // namespace core::solvers::CSRC::symmetric::complex
