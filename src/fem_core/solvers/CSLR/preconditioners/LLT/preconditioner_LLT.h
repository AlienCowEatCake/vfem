#if !defined(SOLVERS_CSLR_PRECONDITIONERS_LLT_H_INCLUDED)
#define SOLVERS_CSLR_PRECONDITIONERS_LLT_H_INCLUDED

#include <cstdio>
#include <cmath>
#include <complex>
#include "../Nothing/preconditioner_Nothing.h"

namespace fem_core { namespace solvers { namespace CSLR { namespace preconditioners {

/**
 * @brief Предобуславливатель LLT, S = L, Q = LT
 */
template<typename val_type, typename ind_type = std::size_t>
class preconditioner_LLT : public preconditioner_Nothing<val_type, ind_type>
{
public:
    preconditioner_LLT(const ind_type * gi, const ind_type * gj, const val_type * di,
                       const val_type * gl, const val_type * gu, ind_type n)
        : preconditioner_Nothing<val_type, ind_type>(gi, gj, di, gl, gu, n)
    {
        /// @todo Сделать LU(sq)
        std::fprintf(stderr, "Non-symmetric LLT is not implemented!\n");
        exit(-1);
    }

    preconditioner_LLT(const ind_type * gi, const ind_type * gj, const val_type * di,
                       const val_type * gg, ind_type n)
        : preconditioner_Nothing<val_type, ind_type>(gi, gj, di, gg, n)
    {
        init();
    }

    virtual std::string get_name() const
    {
        return "LLT";
    }

    virtual void solve_S(const val_type * f, val_type * x) const
    {
        for(ind_type k = 1, k1 = 0; k <= this->m_n; k++, k1++)
        {
            val_type sum = 0.0;

            for(ind_type i = this->m_gi[k1]; i < this->m_gi[k]; i++)
                sum += m_L_gg[i] * x[this->m_gj[i]];

            x[k1] = (f[k1] - sum) / m_L_di[k1];
        }
    }

    virtual void solve_ST(const val_type * f, val_type * x) const
    {
        return solve_S(f, x);
    }

    virtual void solve_Q(const val_type * f, val_type * x) const
    {
//        for(ind_type k = this->m_n, k1 = this->m_n-1; k > 0; k--, k1--)
//        {
//            x[k1] = f[k1] / m_L_di[k1];
//            val_type v_el = x[k1];
//            for(ind_type i = this->m_gi[k1]; i < this->m_gi[k]; i++)
//                const_cast<val_type *>(f)[this->m_gj[i]] -= m_L_gg[i] * v_el;
//        }
        for(ind_type i = 0; i < this->m_n; i++)
            x[i] = f[i];
        for(ind_type k = this->m_n, k1 = this->m_n-1; k > 0; k--, k1--)
        {
            x[k1] /= m_L_di[k1];
            val_type v_el = x[k1];
            for(ind_type i = this->m_gi[k1]; i < this->m_gi[k]; i++)
                x[this->m_gj[i]] -= m_L_gg[i] * v_el;
        }
    }

    virtual void solve_QT(const val_type * f, val_type * x) const
    {
        return solve_Q(f, x);
    }

    virtual void mul_Q(const val_type * f, val_type * x) const
    {
        for(ind_type i = 0; i < this->m_n; i++)
        {
            val_type v_el = f[i];
            x[i] = m_L_di[i] * v_el;
            for(ind_type k = this->m_gi[i], k1 = this->m_gi[i+1]; k < k1; k++)
                x[this->m_gj[k]] += m_L_gg[k] * v_el;
        }
    }

    virtual ~preconditioner_LLT()
    {
        delete [] m_L_di;
        delete [] m_L_gg;
    }

protected:

    val_type * m_L_di, * m_L_gg;

    void init()
    {
        m_L_di = new val_type [this->m_n];
        m_L_gg = new val_type [this->m_gi[this->m_n]];

        for(ind_type i = 0; i < this->m_n; i++)
            m_L_di[i] = this->m_di[i];

        for(ind_type i = 0 ; i < this->m_gi[this->m_n] ; i++)
            m_L_gg[i] = this->m_gl[i];

        make_LLT_decomposition();
    }

    void make_LLT_decomposition()
    {
        val_type sum_d, sum_l;

        for(ind_type k = 0; k < this->m_n; k++)
        {
            sum_d = 0;
            ind_type i_s = this->m_gi[k], i_e = this->m_gi[k+1];

            for(ind_type i = i_s; i < i_e ; i++)
            {
                sum_l = 0;
                ind_type j_s = this->m_gi[this->m_gj[i]], j_e = this->m_gi[this->m_gj[i]+1];

                for(ind_type m = i_s; m < i; m++)
                {
                    for(ind_type j = j_s; j < j_e; j++)
                    {
                        if(this->m_gj[m] == this->m_gj[j])
                        {
                            sum_l += m_L_gg[m] * m_L_gg[j];
                            j_s++;
                        }
                    }
                }
                m_L_gg[i] = (m_L_gg[i] -  sum_l) / m_L_di[this->m_gj[i]];

                sum_d += m_L_gg[i] * m_L_gg[i];
            }
            m_L_di[k] = std::sqrt(m_L_di[k] - sum_d);
        }
    }
};

}}}} // namespace fem_core::solvers::CSLR::preconditioners

#endif // SOLVERS_CSLR_PRECONDITIONERS_LLT_H_INCLUDED
