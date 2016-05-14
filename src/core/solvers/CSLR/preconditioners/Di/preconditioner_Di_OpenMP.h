#if !defined(SOLVERS_CSLR_PRECONDITIONERS_DI_OPENMP_H_INCLUDED)
#define SOLVERS_CSLR_PRECONDITIONERS_DI_OPENMP_H_INCLUDED

#include "../Nothing/preconditioner_Nothing_OpenMP.h"

namespace core { namespace solvers { namespace CSLR { namespace preconditioners {

/**
 * @brief Диагональный предобуславливатель, S = Di, Q = I, распараллеленный через OpenMP
 */
template<typename val_type, typename ind_type = std::size_t>
class preconditioner_Di_OpenMP : public preconditioner_Nothing_OpenMP<val_type, ind_type>
{
    typedef wrappers::omp::omp_int omp_int;

public:
    preconditioner_Di_OpenMP(const ind_type * gi, const ind_type * gj, const val_type * di,
                             const val_type * gl, const val_type * gu, ind_type n)
        : preconditioner_Nothing_OpenMP<val_type, ind_type>(gi, gj, di, gl, gu, n)
    {}

    preconditioner_Di_OpenMP(const ind_type * gi, const ind_type * gj, const val_type * di,
                             const val_type * gg, ind_type n)
        : preconditioner_Nothing_OpenMP<val_type, ind_type>(gi, gj, di, gg, n)
    {}

    virtual std::string get_name() const
    {
        return "Di_OpenMP";
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

}}}} // namespace core::solvers::CSLR::preconditioners

#endif // SOLVERS_CSLR_PRECONDITIONERS_DI_OPENMP_H_INCLUDED
