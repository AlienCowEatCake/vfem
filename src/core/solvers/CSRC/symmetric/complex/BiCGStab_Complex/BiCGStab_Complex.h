#if !defined(SOLVERS_CSRC_SYMMETRIC_COMPLEX_BICGSTAB_H_INCLUDED)
#define SOLVERS_CSRC_SYMMETRIC_COMPLEX_BICGSTAB_H_INCLUDED

#include <cstdlib>
#include <complex>
#include "../../symmetric_solver_interface.h"

namespace core { namespace solvers { namespace CSRC { namespace symmetric { namespace complex {

/**
 * @brief Комплексный BiCGStab
 */
class BiCGStab_Complex : public symmetric_solver_interface<std::complex<double>, std::size_t>
{
public:
    void init(const std::size_t * gi, const std::size_t * gj, const std::complex<double> * di,
              const std::complex<double> * gg, std::size_t n,
              preconditioners::preconditioner_interface<std::complex<double>, std::size_t> * precond);
    void solve(std::complex<double> * solution, const std::complex<double> * rp, double gamma, std::size_t max_iter);

    BiCGStab_Complex();
    ~BiCGStab_Complex();

protected:
    void mul_matrix(const std::complex<double> * f, std::complex<double> * x) const;
    std::complex<double> dot_prod(const std::complex<double> * a, const std::complex<double> * b) const;
    std::complex<double> dot_prod_nocj(const std::complex<double> * a, const std::complex<double> * b) const;
    double dot_prod_self(const std::complex<double> * a) const;

    std::size_t m_n;
    const std::size_t * m_gi, * m_gj;
    const std::complex<double> * m_di, * m_gg;

    std::complex<double> * m_r, * m_r2, * m_v, * m_s, * m_x0, * m_p, * m_t;
};

}}}}} // namespace core::solvers::CSRC::symmetric::complex

#endif // SOLVERS_CSRC_SYMMETRIC_COMPLEX_BICGSTAB_H_INCLUDED
