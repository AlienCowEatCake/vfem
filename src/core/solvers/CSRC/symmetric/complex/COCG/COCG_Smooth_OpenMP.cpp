#if defined(_MSC_VER) && !defined(_CRT_SECURE_NO_WARNINGS)
#define _CRT_SECURE_NO_WARNINGS
#endif
#include "COCG_Smooth_OpenMP.h"
#include <cstdio>
#include <cmath>
#include "../../../../../utils/fpu.h"
#include "../../../../../wrappers/omp_wrapper.h"

namespace core { namespace solvers { namespace CSRC { namespace symmetric { namespace complex {

typedef wrappers::omp::omp_int omp_int;

void COCG_Smooth_OpenMP::init(const std::size_t * gi, const std::size_t * gj, const std::complex<double> * di,
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
    delete [] m_mv_tmp;
    delete [] m_mv_ind;

    m_r = new std::complex<double> [m_n];
    m_z = new std::complex<double> [m_n];
    m_p = new std::complex<double> [m_n];
    m_s = new std::complex<double> [m_n];
    m_xs = new std::complex<double> [m_n];
    m_rs = new std::complex<double> [m_n];
    m_mv_tmp = new std::complex<double> [n * static_cast<std::size_t>(m_num_threads - 1)];

    // Индексы начала/конца в f и x
    m_mv_ind = new std::size_t [m_num_threads + 1];
    // Начинает первый в начале
    m_mv_ind[0] = 0;
    // Заканчивает последний в конце
    m_mv_ind[m_num_threads] = n;
    // Выравнивание по строкам
    for(std::size_t i = 0, curr_ind = 1, curr_size = gi[n] / static_cast<std::size_t>(m_num_threads); i < n; i++)
    {
        // Если где-то между текущей и следующей строкой проходит разбиение
        if(gi[i + 1] >= curr_size && gi[i] <= curr_size)
        {
            // Найдем куда ближе сдвигать - к текущей или к следующей
            // Если ближе к следующей - сделаем ее текущей
            if(gi[i + 1] - curr_size < curr_size - gi[i])
                i++;
            // Сдвигаемся
            m_mv_ind[curr_ind] = i;
            curr_size += gi[n] / static_cast<std::size_t>(m_num_threads);
            curr_ind++;
        }
    }
}

std::complex<double> COCG_Smooth_OpenMP::dot_prod_nocj(const std::complex<double> * a, const std::complex<double> * b) const
{
    double d_p_real = 0.0;
    double d_p_imag = 0.0;
#pragma omp parallel for reduction(+ : d_p_real, d_p_imag)
    for(omp_int i = 0; i < static_cast<omp_int>(m_n); i++)
    {
        std::complex<double> c = a[i] * b[i];
        d_p_real += c.real();
        d_p_imag += c.imag();
    }
    return std::complex<double>(d_p_real, d_p_imag);
}

double COCG_Smooth_OpenMP::dot_prod_self(const std::complex<double> * a) const
{
    double d_p = 0.0;
#pragma omp parallel for reduction(+ : d_p)
    for(omp_int i = 0; i < static_cast<omp_int>(m_n); i++)
    {
        double re = a[i].real();
        double im = a[i].imag();
        d_p += re * re + im * im;
    }
    return d_p;
}

double COCG_Smooth_OpenMP::dot_prod_real(const std::complex<double> * a, const std::complex<double> * b) const
{
    double d_p = 0.0;
#pragma omp parallel for reduction(+ : d_p)
    for(omp_int i = 0; i < static_cast<omp_int>(m_n); i++)
        d_p += a[i].real() * b[i].real() + a[i].imag() * b[i].imag();
    return d_p;
}

void COCG_Smooth_OpenMP::mul_matrix(const std::complex<double> * f, std::complex<double> * x) const
{
    std::size_t nt1 = static_cast<std::size_t>(m_num_threads - 1);

#pragma omp parallel num_threads(m_num_threads)
    {
#pragma omp for nowait
        for(omp_int i = 0; i < static_cast<omp_int>(m_n); i++)
            x[i] = std::complex<double>(0.0, 0.0);

#pragma omp for
        for(omp_int j = 0; j < static_cast<omp_int>(m_n * nt1); j++)
            m_mv_tmp[j] = std::complex<double>(0.0, 0.0);

        omp_int rank = static_cast<omp_int>(omp_get_thread_num());
        omp_int i_beg = static_cast<omp_int>(m_mv_ind[rank]);
        omp_int i_end = static_cast<omp_int>(m_mv_ind[rank + 1]);
        for(omp_int i = i_beg; i < i_end; i++)
        {
            if(rank == 0)
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
            else
            {
                std::size_t adr = (rank - 1) * m_n;
                std::size_t adr_i = adr + i;
                std::complex<double> v_el = f[i];
                m_mv_tmp[adr_i] = m_di[i] * v_el;
                for(std::size_t k = m_gi[i], k1 = m_gi[i + 1]; k < k1; k++)
                {
                    std::size_t j = m_gj[k];
                    m_mv_tmp[adr_i] += m_gg[k] * f[j];
                    m_mv_tmp[adr + j] += m_gg[k] * v_el;
                }
            }
        }
#pragma omp barrier

#pragma omp for
        for(omp_int i = 0; i < static_cast<omp_int>(m_n); i++)
        {
            for(std::size_t j = 0; j < nt1; j++)
                x[i] += m_mv_tmp[j * m_n + i];
        }
    }
}

void COCG_Smooth_OpenMP::solve_SQ(const std::complex<double> * f, std::complex<double> * x) const
{
    m_precond->solve_S(f, x);
    m_precond->solve_Q(x, x);
}

void COCG_Smooth_OpenMP::solve(std::complex<double> * solution, const std::complex<double> * rp,
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
        m_xs[i] = m_x0[i];
        m_rs[i] = m_r[i] = m_rp[i] - m_r[i];
    }

    solve_SQ(m_r, m_z);
#pragma omp parallel for
    for(omp_int i = 0; i < static_cast<omp_int>(m_n); i++)
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
#pragma omp parallel for
            for(omp_int i = 0; i < static_cast<omp_int>(m_n); i++)
                solution[i] = m_xs[i];
            return;
        }

        double residual = discr / rp_norm;
//        if(iter%10 == 0)
        {
            printf("COCG_Smooth_OpenMP<%s> [%d] Residual:\t%5lu\t%.3e\r", precond_name, m_num_threads, static_cast<unsigned long>(iter), sqrt(residual));
            fflush(stdout);
        }

        if(residual > eps)
        {
            mul_matrix(m_z, m_s);

            alpha = prod_1 / dot_prod_nocj(m_s, m_z);

#pragma omp parallel for
            for(omp_int i = 0; i < static_cast<omp_int>(m_n); i++)
            {
                m_x0[i] += alpha * m_z[i];
                m_r[i] -= alpha * m_s[i];
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

#pragma omp parallel for
            for(omp_int i = 0; i < static_cast<omp_int>(m_n); i++)
                m_z[i] = m_p[i] + beta * m_z[i];
        }
        else
            finished = true;
    }

//    mul_matrix(m_xs, m_r);
//#pragma omp parallel for
//    for(omp_int i = 0; i < static_cast<omp_int>(m_n); i++)
//        m_r[i] = m_rp[i] - m_r[i];
//    discr = dot_prod_self(m_r);
    printf("COCG_Smooth_OpenMP<%s> [%d] Residual:\t%5lu\t%.3e\n", precond_name_str.c_str(), m_num_threads, static_cast<unsigned long>(iter) - 1, sqrt(discr / rp_norm));

    if(iter >= max_iter)
        printf("Soulution can`t found, iteration limit exceeded!\n");

#pragma omp parallel for
    for(omp_int i = 0; i < static_cast<omp_int>(m_n); i++)
        solution[i] = m_xs[i];
}

COCG_Smooth_OpenMP::COCG_Smooth_OpenMP()
{
    wrappers::omp::wrapper_omp_set_env_max_threads();
    m_num_threads = omp_get_max_threads();
    m_r = m_x0 = m_z = m_p = m_s = m_xs = m_rs = NULL;
    m_mv_tmp = NULL;
    m_mv_ind = NULL;
}

COCG_Smooth_OpenMP::~COCG_Smooth_OpenMP()
{
    delete [] m_r;
    delete [] m_z;
    delete [] m_p;
    delete [] m_s;
    delete [] m_xs;
    delete [] m_rs;
    delete [] m_mv_tmp;
    delete [] m_mv_ind;
}

}}}}} // namespace core::solvers::CSRC::symmetric::complex
