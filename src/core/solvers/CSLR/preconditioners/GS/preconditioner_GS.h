#if !defined(SOLVERS_CSLR_PRECONDITIONERS_GS_H_INCLUDED)
#define SOLVERS_CSLR_PRECONDITIONERS_GS_H_INCLUDED

#include "../Nothing/preconditioner_Nothing.h"

namespace core { namespace solvers { namespace CSLR { namespace preconditioners {

/**
 * @brief Предобуславливатель типа Гаусса-Зейделя, S = (D + L)^-1 * D * (D + U)^-1, Q = I
 */
template<typename val_type, typename ind_type = std::size_t>
class preconditioner_GS : public preconditioner_Nothing<val_type, ind_type>
{
public:
    preconditioner_GS(const ind_type * gi, const ind_type * gj, const val_type * di,
                      const val_type * gl, const val_type * gu, ind_type n)
        : preconditioner_Nothing<val_type, ind_type>(gi, gj, di, gl, gu, n)
    {
        init();
    }

    preconditioner_GS(const ind_type * gi, const ind_type * gj, const val_type * di,
                      const val_type * gg, ind_type n)
        : preconditioner_Nothing<val_type, ind_type>(gi, gj, di, gg, n)
    {
        init();
    }

    virtual std::string get_name() const
    {
        return "GS";
    }

    virtual void solve_S(const val_type * f, val_type * x) const
    {
        val_type w = 1;

        for(ind_type i = 0; i < this->m_n; i++)
            this->m_tmp[i] = f[i];

        // (D + L)^-1 tmp = x
        for(ind_type k = 1, k1 = 0; k <= this->m_n; k++, k1++)
        {
            val_type sum = 0;
            for(ind_type i = this->m_gi[k1]; i < this->m_gi[k]; i++)
                sum += this->m_gl[i] * x[this->m_gj[i]] * w;
            x[k1] = (this->m_tmp[k1] - sum) / this->m_di[k1];
        }

        // tmp = D * x
        for(ind_type i = 0; i < this->m_n; i++)
            this->m_tmp[i] = x[i] * this->m_di[i];

        // (D + U)^-1 * tmp = x
        for(size_t k = this->m_n, k1 = this->m_n - 1; k > 0; k--, k1--)
        {
            x[k1] = this->m_tmp[k1] / this->m_di[k1];
            val_type v_el = x[k1];
            for(ind_type i = this->m_gi[k1]; i < this->m_gi[k]; i++)
                this->m_tmp[this->m_gj[i]] -= this->m_gu[i] * v_el * w;
        }
    }

    virtual ~preconditioner_GS()
    {
        delete [] m_tmp;
    }

protected:

    mutable val_type * m_tmp;

    void init()
    {
        m_tmp = new val_type [this->m_n];
    }
};

}}}} // namespace core::solvers::CSLR::preconditioners

#endif // SOLVERS_CSLR_PRECONDITIONERS_GS_H_INCLUDED
