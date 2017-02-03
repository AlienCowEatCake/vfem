#if !defined(SOLVERS_CSLR_SYMMETRIC_COMPLEX_COCG_SMOOTH_OPENMP_H_INCLUDED)
#define SOLVERS_CSLR_SYMMETRIC_COMPLEX_COCG_SMOOTH_OPENMP_H_INCLUDED

#include <cstdlib>
#include <complex>
#include "../../symmetric_solver_interface.h"

namespace fem_core { namespace solvers { namespace CSLR { namespace symmetric { namespace complex {

/**
 * @brief COCG со сглаживанием невязки, распараллеленный через OpenMP
 */
class COCG_Smooth_OpenMP : public symmetric_solver_interface<std::complex<double>, std::size_t>
{
public:
    void init(const std::size_t * gi, const std::size_t * gj, const std::complex<double> * di,
              const std::complex<double> * gg, std::size_t n,
              preconditioners::preconditioner_interface<std::complex<double>, std::size_t> * precond);
    void solve(std::complex<double> * solution, const std::complex<double> * rp,
               double eps, size_t max_iter);

    COCG_Smooth_OpenMP();
    ~COCG_Smooth_OpenMP();

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
    std::complex<double> * m_mv_tmp;
    std::size_t * m_mv_ind;
};

}}}}} // namespace fem_core::solvers::CSLR::symmetric::complex

#endif // SOLVERS_CSLR_SYMMETRIC_COMPLEX_COCG_SMOOTH_OPENMP_H_INCLUDED
