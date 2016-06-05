#if !defined(SOLVERS_CSLR_SYMMETRIC_COMPLEX_COCG_SMOOTH_H_INCLUDED)
#define SOLVERS_CSLR_SYMMETRIC_COMPLEX_COCG_SMOOTH_H_INCLUDED

#include "COCG.h"

namespace core { namespace solvers { namespace CSLR { namespace symmetric { namespace complex {

/**
 * @brief COCG со сглаживанием невязки
 */
class COCG_Smooth : public COCG
{
public:
    void init(const std::size_t * gi, const std::size_t * gj, const std::complex<double> * di,
              const std::complex<double> * gg, std::size_t n,
              preconditioners::preconditioner_interface<std::complex<double>, std::size_t> * precond);
    void solve(std::complex<double> * solution, const std::complex<double> * rp,
               double eps, size_t max_iter);

    COCG_Smooth();
    ~COCG_Smooth();

protected:
    double dot_prod_real(const std::complex<double> * a, const std::complex<double> * b) const;

    std::complex<double> * m_xs, * m_rs;
};

}}}}} // namespace core::solvers::CSLR::symmetric::complex

#endif // SOLVERS_CSLR_SYMMETRIC_COMPLEX_COCG_SMOOTH_H_INCLUDED
