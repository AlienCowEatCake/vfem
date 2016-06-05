#if !defined(SOLVERS_CSLR_SYMMETRIC_COMPLEX_COCG_SMOOTH_MKL_H_INCLUDED)
#define SOLVERS_CSLR_SYMMETRIC_COMPLEX_COCG_SMOOTH_MKL_H_INCLUDED

#include <cstdlib>
#include <complex>
#include "../../symmetric_solver_interface.h"
#include "../../../../../wrappers/mkl_wrapper.h"

namespace core { namespace solvers { namespace CSLR { namespace symmetric { namespace complex {

/**
 * @brief COCG со сглаживанием невязки, использующий библиотеку MKL
 */
class COCG_Smooth_MKL : public symmetric_solver_interface<std::complex<double>, std::size_t>
{
public:
    void init(const std::size_t * gi, const std::size_t * gj, const std::complex<double> * di,
              const std::complex<double> * gg, std::size_t n,
              preconditioners::preconditioner_interface<std::complex<double>, std::size_t> * precond);
    void solve(std::complex<double> * solution, const std::complex<double> * rp,
               double eps, size_t max_iter);

    COCG_Smooth_MKL();
    ~COCG_Smooth_MKL();

protected:
    void mul_matrix(const std::complex<double> * f, std::complex<double> * x) const;
    void solve_SQ(const std::complex<double> * f, std::complex<double> * x) const;
    std::complex<double> dot_prod_nocj(const std::complex<double> * a, const std::complex<double> * b) const;
    double dot_prod_self(const std::complex<double> * a) const;
    double dot_prod_real(const std::complex<double> * a, const std::complex<double> * b) const;

    std::size_t m_n;
    const std::size_t * m_gi, * m_gj;
    const std::complex<double> * m_di, * m_gg, * m_rp;
    std::complex<double> * m_r, * m_x0, * m_z, * m_p, * m_s, * m_xs, * m_rs;
    preconditioners::preconditioner_interface<std::complex<double>, std::size_t> * m_precond;

    int m_num_threads;
    MKL_INT m_m;
    MKL_INT * m_ia, * m_ja;
    MKL_Complex16 * m_aa;
};

}}}}} // namespace core::solvers::CSLR::symmetric::complex

#endif // SOLVERS_CSLR_SYMMETRIC_COMPLEX_COCG_SMOOTH_MKL_H_INCLUDED
