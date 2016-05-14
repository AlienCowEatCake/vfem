#if !defined(SOLVERS_CSLR_SYMMETRIC_COMPLEX_BICGSTAB_SMOOTH_H_INCLUDED)
#define SOLVERS_CSLR_SYMMETRIC_COMPLEX_BICGSTAB_SMOOTH_H_INCLUDED

#include "BiCGStab_Complex.h"

namespace core { namespace solvers { namespace CSLR { namespace symmetric { namespace complex {

/**
 * @brief Комплексный BiCGStab со сглаживанием невязки
 */
class BiCGStab_Complex_Smooth : public BiCGStab_Complex
{
public:
    void init(const std::size_t * gi, const std::size_t * gj, const std::complex<double> * di,
              const std::complex<double> * gg, std::size_t n,
              preconditioners::preconditioner_interface<std::complex<double>, std::size_t> * precond);
    void solve(std::complex<double> * solution, const std::complex<double> * rp, double gamma, std::size_t max_iter);

    BiCGStab_Complex_Smooth();
    ~BiCGStab_Complex_Smooth();

protected:
    double dot_prod_real(const std::complex<double> * a, const std::complex<double> * b) const;

    std::complex<double> * m_xs, * m_rs;
};

}}}}} // namespace core::solvers::CSLR::symmetric::complex

#endif // SOLVERS_CSLR_SYMMETRIC_COMPLEX_BICGSTAB_SMOOTH_H_INCLUDED
