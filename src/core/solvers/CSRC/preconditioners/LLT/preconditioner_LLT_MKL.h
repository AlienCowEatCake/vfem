#if !defined(SOLVERS_CSRC_PRECONDITIONERS_LLT_MKL_H_INCLUDED)
#define SOLVERS_CSRC_PRECONDITIONERS_LLT_MKL_H_INCLUDED

#include <cstdio>
#include <cmath>
#include <complex>
#include "../Nothing/preconditioner_Nothing_MKL.h"
#include "../../../../wrappers/omp_wrapper.h"
#include "../../../../wrappers/mkl_wrapper.h"

namespace core { namespace solvers { namespace CSRC { namespace preconditioners {

/**
 * @brief Предобуславливатель LLT, S = L, Q = LT, использующий библиотеку MKL
 */
template<typename ind_type = std::size_t>
class preconditioner_LLT_MKL : public preconditioner_Nothing_MKL<std::complex<double>, ind_type>
{
    typedef wrappers::omp::omp_int omp_int;

public:
    preconditioner_LLT_MKL(const ind_type * gi, const ind_type * gj, const std::complex<double> * di,
                           const std::complex<double> * gl, const std::complex<double> * gu, ind_type n)
        : preconditioner_Nothing_MKL<std::complex<double>, ind_type>(gi, gj, di, gl, gu, n)
    {
        /// @todo Сделать LU(sq)
        std::fprintf(stderr, "Non-symmetric LLT_MKL is not implemented!\n");
        exit(-1);
    }

    preconditioner_LLT_MKL(const ind_type * gi, const ind_type * gj, const std::complex<double> * di,
                           const std::complex<double> * gg, ind_type n)
        : preconditioner_Nothing_MKL<std::complex<double>, ind_type>(gi, gj, di, gg, n)
    {
        init();
    }

    virtual std::string get_name() const
    {
        return "LLT_MKL";
    }

    virtual void solve_S(const std::complex<double> * f, std::complex<double> * x) const
    {
        const MKL_Complex16 * ff = reinterpret_cast<const MKL_Complex16 *>(f);
        MKL_Complex16 * xx = reinterpret_cast<MKL_Complex16 *>(x);
        mkl_zcsrtrsv("L", "N", "N", &m_mkl, L_aa, ia, ja, ff, xx);
    }

    virtual void solve_ST(const std::complex<double> * f, std::complex<double> * x) const
    {
        return solve_S(f, x);
    }

    virtual void solve_Q(const std::complex<double> * f, std::complex<double> * x) const
    {
        const MKL_Complex16 * ff = reinterpret_cast<const MKL_Complex16 *>(f);
        MKL_Complex16 * xx = reinterpret_cast<MKL_Complex16 *>(x);
        mkl_zcsrtrsv("L", "T", "N", &m_mkl, L_aa, ia, ja, ff, xx);
    }

    virtual void solve_QT(const std::complex<double> * f, std::complex<double> * x) const
    {
        return solve_Q(f, x);
    }

    virtual void mul_Q(const std::complex<double> * f, std::complex<double> * x) const
    {
        MKL_Complex16 * in_v = reinterpret_cast<MKL_Complex16 *>(const_cast<std::complex<double> *>(f));
        MKL_Complex16 * out_v = reinterpret_cast<MKL_Complex16 *>(const_cast<std::complex<double> *>(x));
        mkl_zcsrsymv("L", &m_mkl, L_aa, ia, ja, in_v, out_v);
    }

    virtual ~preconditioner_LLT_MKL()
    {
        delete [] L_aa;
        delete [] ia;
        delete [] ja;
    }

protected:

    MKL_Complex16 * L_aa;
    MKL_INT m_mkl;
    MKL_INT * ia, * ja;

    void init()
    {
        const std::complex<double> * gg = this->m_gl;
        m_mkl = static_cast<MKL_INT>(this->m_n);
        ia = new MKL_INT[this->m_n + 1];
        ia[0] = 1;
        for(MKL_INT i = 1; i <= m_mkl; i++)
            ia[i] = ia[i - 1] + static_cast<MKL_INT>(this->m_gi[i] - this->m_gi[i - 1]) + 1;
        ja = new MKL_INT[ia[this->m_n]];
        L_aa = new MKL_Complex16[ia[this->m_n]];
#pragma omp parallel for
        for(omp_int i = 0; i < static_cast<omp_int>(this->m_n); i++)
        {
            ind_type i_b = this->m_gi[i], i_e = this->m_gi[i + 1];
            for(ind_type j = i_b; j < i_e; j++)
            {
                ind_type ind = ia[i] + (j - this->m_gi[i]) - 1;
                ja[ind] = static_cast<MKL_INT>(this->m_gj[j]) + 1;
                L_aa[ind].real = gg[j].real();
                L_aa[ind].imag = gg[j].imag();
            }
            ind_type ind = ia[i + 1] - 2;
            ja[ind] = i + 1;
            L_aa[ind].real = this->m_di[i].real();
            L_aa[ind].imag = this->m_di[i].imag();
        }
        make_LLT_MKL_decomposition();
    }

    void make_LLT_MKL_decomposition()
    {
        std::complex<double> * L_gg = reinterpret_cast<std::complex<double> *>(L_aa);

        std::complex<double> sum_d, sum_l;

        for(ind_type k = 0; k < this->m_n; k++)
        {
            sum_d = 0;
            ind_type i_s = ia[k] - 1, i_e = ia[k+1] - 2;
            for(ind_type i = i_s; i < i_e ; i++)
            {
                sum_l = 0;
                ind_type j_s = ia[ja[i]-1]-1, j_e = ia[ja[i]]-2;

                for(ind_type m = i_s; m < i; m++)
                {
                    for(ind_type j = j_s; j < j_e; j++)
                    {
                        if(ja[m] == ja[j])
                        {
                            sum_l += L_gg[m] * L_gg[j];
                            j_s++;
                        }
                    }
                }
                L_gg[i] = (L_gg[i] -  sum_l) / L_gg[j_e];

                sum_d += L_gg[i] * L_gg[i];
            }
            ind_type i_di = ia[k+1]-2;
            L_gg[i_di] = sqrt(L_gg[i_di] - sum_d);
        }
    }
};

}}}} // namespace core::solvers::CSRC::preconditioners

#endif // SOLVERS_CSRC_PRECONDITIONERS_LLT_MKL_H_INCLUDED
