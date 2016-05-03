#if defined(_MSC_VER) && !defined(_CRT_SECURE_NO_WARNINGS)
#define _CRT_SECURE_NO_WARNINGS
#endif
#include "GMRES_Complex_OpenMP.h"
#include <cstdio>
#include <cmath>
#include "../../../../../utils/fpu.h"
#include "../../../../../wrappers/omp_wrapper.h"

namespace core { namespace solvers { namespace CSRC { namespace symmetric { namespace complex {

typedef wrappers::omp::omp_int omp_int;

void GMRES_Complex_OpenMP::init(const std::size_t * gi, const std::size_t * gj, const std::complex<double> * di,
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
            m_m = m_m_curr = (std::size_t)temp_m;
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

    delete [] m_mv_tmp;
    delete [] m_mv_ind;
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

std::complex<double> GMRES_Complex_OpenMP::dot_prod(const std::complex<double> * a, const std::complex<double> * b) const
{
    double d_p_real = 0.0;
    double d_p_imag = 0.0;
#pragma omp parallel for reduction(+ : d_p_real, d_p_imag)
    for(omp_int i = 0; i < static_cast<omp_int>(m_n); i++)
    {
        std::complex<double> c = conj(a[i]) * b[i];
        d_p_real += c.real();
        d_p_imag += c.imag();
    }
    return std::complex<double>(d_p_real, d_p_imag);
}

double GMRES_Complex_OpenMP::dot_prod_self(const std::complex<double> * a) const
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

void GMRES_Complex_OpenMP::mul_matrix(const std::complex<double> * f, std::complex<double> * x) const
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

void GMRES_Complex_OpenMP::solve_QAS(const std::complex<double> * f, std::complex<double> * x, std::complex<double> * tmp) const
{
    copy_vec(f, x);
    m_precond->solve_Q(x, x);
    mul_matrix(x, tmp);
    m_precond->solve_S(tmp, x);
}

void GMRES_Complex_OpenMP::copy_vec(const std::complex<double> * f, std::complex<double> * x) const
{
#pragma omp parallel for
    for(omp_int i = 0; i < static_cast<omp_int>(m_n); i++)
        x[i] = f[i];
}

void GMRES_Complex_OpenMP::solve(std::complex<double> * solution, const std::complex<double> * rp,
                                 double eps, std::size_t max_iter)
{
    m_rp = rp;
    std::string precond_name_str = m_precond->get_name();
    const char * precond_name = precond_name_str.c_str();

    m_x0 = new std::complex<double> [m_n];
    m_precond->mul_Q(solution, m_x0);

    mul_matrix(solution, m_r);
#pragma omp parallel for
    for(omp_int i = 0; i < static_cast<omp_int>(m_n); i++)
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
            printf("GMRES_Complex_OpenMP<%s>(%lu) [%d] Residual:\t%5lu\t%.3e\r", precond_name, (unsigned long)m_m, m_num_threads, (unsigned long)iter, residual);
            fflush(stdout);
        }

        if(residual > eps && fabs(residual - residual_prev) / (residual) > 1e-15)
        {
            // V
            #pragma omp parallel for
            for(omp_int i = 0; i < static_cast<omp_int>(m_n); i++)
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
                {
                    #pragma omp parallel for
                    for(omp_int i = 0; i < static_cast<omp_int>(m_n); i++)
                        m_w[i] -= m_VT[lambda][i] * m_H[lambda][mu];
                }

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
                {
                    #pragma omp parallel for
                    for(omp_int i = 0; i < static_cast<omp_int>(m_n); i++)
                        m_VT[mu + 1][i] = m_w[i] / vnorm;
                }
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
            #pragma omp parallel for
            for(omp_int i = 0; i < static_cast<omp_int>(m_n); i++)
                for(std::size_t j = 0; j < m_m_curr; j++)
                    m_x0[i] += m_VT[j][i] * m_d[j];

            // r
            copy_vec(m_x0, m_t);
            m_precond->solve_Q(m_t, m_t);
            mul_matrix(m_t, m_r);
            #pragma omp parallel for
            for(omp_int i = 0; i < static_cast<omp_int>(m_n); i++)
                m_t[i] = m_rp[i] - m_r[i];
            m_precond->solve_S(m_t, m_r);
        }
        else
            finished = true;
    }

    m_precond->solve_Q(m_x0, solution);
    delete [] m_x0;

//    mul_matrix(solution, m_r);
//#pragma omp parallel for
//    for(omp_int i = 0; i < static_cast<omp_int>(m_n); i++)
//        m_r[i] = m_rp[i] - m_r[i];
//    discr = sqrt(dot_prod_self(m_r));
//    rp_norm = sqrt(dot_prod_self(m_rp));
    printf("GMRES_Complex_OpenMP<%s>(%lu) [%d] Residual:\t%5lu\t%.3e\n", precond_name_str.c_str(), (unsigned long)m_m, m_num_threads, (unsigned long)iter - 1, discr / rp_norm);

    if(iter >= max_iter)
        printf("Soulution can`t found, iteration limit exceeded!\n");
    else if(residual > eps)
        printf("Soulution can`t found, stagnation detected!\n");
}

GMRES_Complex_OpenMP::GMRES_Complex_OpenMP()
{
    m_d = m_r = m_w = m_t = m_x0 = NULL;
    m_H = m_VT = NULL;
    m_m = m_m_curr = 0;

    wrappers::omp::wrapper_omp_set_env_max_threads();
    m_num_threads = omp_get_max_threads();
    m_mv_tmp = NULL;
    m_mv_ind = NULL;
}

GMRES_Complex_OpenMP::~GMRES_Complex_OpenMP()
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

    delete [] m_mv_tmp;
    delete [] m_mv_ind;
}

}}}}} // namespace core::solvers::CSRC::symmetric::complex
