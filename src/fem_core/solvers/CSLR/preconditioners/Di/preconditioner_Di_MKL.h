#if !defined(SOLVERS_CSLR_PRECONDITIONERS_DI_MKL_H_INCLUDED)
#define SOLVERS_CSLR_PRECONDITIONERS_DI_MKL_H_INCLUDED

#include "../Nothing/preconditioner_Nothing_MKL.h"

namespace fem_core { namespace solvers { namespace CSLR { namespace preconditioners {

/**
 * @brief Диагональный предобуславливатель, S = Di, Q = I, использующий библиотеку MKL
 */
template<typename val_type, typename ind_type = std::size_t>
class preconditioner_Di_MKL : public preconditioner_Nothing_MKL<val_type, ind_type>
{
    typedef wrappers::omp::omp_int omp_int;

public:
    preconditioner_Di_MKL(const ind_type * gi, const ind_type * gj, const val_type * di,
                          const val_type * gl, const val_type * gu, ind_type n)
        : preconditioner_Nothing_MKL<val_type, ind_type>(gi, gj, di, gl, gu, n)
    {}

    preconditioner_Di_MKL(const ind_type * gi, const ind_type * gj, const val_type * di,
                          const val_type * gg, ind_type n)
        : preconditioner_Nothing_MKL<val_type, ind_type>(gi, gj, di, gg, n)
    {}

    virtual std::string get_name() const
    {
        return "Di_MKL";
    }

    virtual void solve_S(const val_type * f, val_type * x) const
    {
#pragma omp parallel for
        for(omp_int k = 0; k < static_cast<omp_int>(this->m_n); k++)
            x[k] = f[k] / this->m_di[k];
    }

    virtual void solve_ST(const val_type * f, val_type * x) const
    {
#pragma omp parallel for
        for(omp_int k = 0; k < static_cast<omp_int>(this->m_n); k++)
            x[k] = f[k] / this->m_di[k];
    }
};

}}}} // namespace fem_core::solvers::CSLR::preconditioners

#endif // SOLVERS_CSLR_PRECONDITIONERS_DI_MKL_H_INCLUDED
