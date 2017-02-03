#if !defined(SOLVERS_CSLR_SYMMETRIC_COMPLEX_COCR_H_INCLUDED)
#define SOLVERS_CSLR_SYMMETRIC_COMPLEX_COCR_H_INCLUDED

#include <cstdlib>
#include <complex>
#include "../../symmetric_solver_interface.h"

namespace fem_core { namespace solvers { namespace CSLR { namespace symmetric { namespace complex {

/**
 * @brief COCR
 */
class COCR : public symmetric_solver_interface<std::complex<double>, std::size_t>
{
public:
    void init(const std::size_t * gi, const std::size_t * gj, const std::complex<double> * di,
              const std::complex<double> * gg, std::size_t n,
              preconditioners::preconditioner_interface<std::complex<double>, std::size_t> * precond);
    void solve(std::complex<double> * solution, const std::complex<double> * rp,
               double eps, std::size_t max_iter);

    COCR();
    ~COCR();

protected:
    void mul_matrix(const std::complex<double> * f, std::complex<double> * x) const;
    void solve_SQ(const std::complex<double> * f, std::complex<double> * x) const;
    std::complex<double> dot_prod_nocj(const std::complex<double> * a, const std::complex<double> * b) const;
    double dot_prod_self(const std::complex<double> * a) const;

    std::size_t m_n;
    const std::size_t * m_gi, * m_gj;
    const std::complex<double> * m_di, * m_gg, * m_rp;
    std::complex<double> * m_r, * m_x0, * m_z, * m_p, * m_s, * m_w, * m_a;
    preconditioners::preconditioner_interface<std::complex<double>, std::size_t> * m_precond;
};

}}}}} // namespace fem_core::solvers::CSLR::symmetric::complex

#endif // SOLVERS_CSLR_SYMMETRIC_COMPLEX_COCR_H_INCLUDED
