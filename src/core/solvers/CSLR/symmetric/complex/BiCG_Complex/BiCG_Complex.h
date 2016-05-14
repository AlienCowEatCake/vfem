#if !defined(SOLVERS_CSLR_SYMMETRIC_COMPLEX_BICG_H_INCLUDED)
#define SOLVERS_CSLR_SYMMETRIC_COMPLEX_BICG_H_INCLUDED

#include <cstdlib>
#include <complex>
#include "../../symmetric_solver_interface.h"

namespace core { namespace solvers { namespace CSLR { namespace symmetric { namespace complex {

/**
 * @brief Комплексный BiCG
 */
class BiCG_Complex : public symmetric_solver_interface<std::complex<double>, std::size_t>
{
public:
    void init(const std::size_t * gi, const std::size_t * gj, const std::complex<double> * di,
              const std::complex<double> * gg, std::size_t n,
              preconditioners::preconditioner_interface<std::complex<double>, std::size_t> * precond);
    void solve(std::complex<double> * solution, const std::complex<double> * rp, double gamma, std::size_t max_iter);

    BiCG_Complex();
    ~BiCG_Complex();

protected:
    void mul_matrix(const std::complex<double> * f, std::complex<double> * x) const;
    std::complex<double> dot_prod(const std::complex<double> * a, const std::complex<double> * b) const;
    std::complex<double> dot_prod_nocj(const std::complex<double> * a, const std::complex<double> * b) const;
    double dot_prod_self(const std::complex<double> * a) const;

    std::size_t m_n;
    const std::size_t * m_gi, * m_gj;
    const std::complex<double> * m_di, * m_gg;

    std::complex<double> * m_r, * m_z, * m_s, * m_x0, * m_p, * m_t, * m_t1;
};

}}}}} // namespace core::solvers::CSLR::symmetric::complex

#endif // SOLVERS_CSLR_SYMMETRIC_COMPLEX_BICG_H_INCLUDED
