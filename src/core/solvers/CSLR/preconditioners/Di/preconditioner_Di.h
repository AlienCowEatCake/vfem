#if !defined(SOLVERS_CSLR_PRECONDITIONERS_DI_H_INCLUDED)
#define SOLVERS_CSLR_PRECONDITIONERS_DI_H_INCLUDED

#include "../Nothing/preconditioner_Nothing.h"

namespace core { namespace solvers { namespace CSLR { namespace preconditioners {

/**
 * @brief Диагональный предобуславливатель, S = Di, Q = I
 */
template<typename val_type, typename ind_type = std::size_t>
class preconditioner_Di : public preconditioner_Nothing<val_type, ind_type>
{
public:
    preconditioner_Di(const ind_type * gi, const ind_type * gj, const val_type * di,
                      const val_type * gl, const val_type * gu, ind_type n)
        : preconditioner_Nothing<val_type, ind_type>(gi, gj, di, gl, gu, n)
    {}

    preconditioner_Di(const ind_type * gi, const ind_type * gj, const val_type * di,
                      const val_type * gg, ind_type n)
        : preconditioner_Nothing<val_type, ind_type>(gi, gj, di, gg, n)
    {}

    virtual std::string get_name() const
    {
        return "Di";
    }

    virtual void solve_S(const val_type * f, val_type * x) const
    {
        for(ind_type k = 0; k < this->m_n; k++)
            x[k] = f[k] / this->m_di[k];
    }

    virtual void solve_ST(const val_type * f, val_type * x) const
    {
        for(ind_type k = 0; k < this->m_n; k++)
            x[k] = f[k] / this->m_di[k];
    }
};

}}}} // namespace core::solvers::CSLR::preconditioners

#endif // SOLVERS_CSLR_PRECONDITIONERS_DI_H_INCLUDED
