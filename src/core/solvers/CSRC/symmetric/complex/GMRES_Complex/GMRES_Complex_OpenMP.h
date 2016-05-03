#if !defined(SOLVERS_CSRC_SYMMETRIC_COMPLEX_GMRES_OPENMP_H_INCLUDED)
#define SOLVERS_CSRC_SYMMETRIC_COMPLEX_GMRES_OPENMP_H_INCLUDED

#include <cstdlib>
#include <complex>
#include "../../symmetric_solver_interface.h"

namespace core { namespace solvers { namespace CSRC { namespace symmetric { namespace complex {

/**
 * @brief Комплексный GMRES, распараллеленный через OpenMP
 */
class GMRES_Complex_OpenMP : public symmetric_solver_interface<std::complex<double>, std::size_t>
{
public:
    void init(const std::size_t * gi, const std::size_t * gj, const std::complex<double> * di,
              const std::complex<double> * gg, std::size_t n,
              preconditioners::preconditioner_interface<std::complex<double>, std::size_t> * precond);
    void solve(std::complex<double> * solution, const std::complex<double> * rp,
               double eps, std::size_t max_iter);

    GMRES_Complex_OpenMP();
    ~GMRES_Complex_OpenMP();

protected:
    void mul_matrix(const std::complex<double> * f, std::complex<double> * x) const;
    void solve_QAS(const std::complex<double> * f, std::complex<double> * x, std::complex<double> * tmp) const;
    std::complex<double> dot_prod(const std::complex<double> * a, const std::complex<double> * b) const;
    double dot_prod_self(const std::complex<double> * a) const;
    void copy_vec(const std::complex<double> * f, std::complex<double> * x) const;

    std::size_t m_n;
    const std::size_t * m_gi, * m_gj;
    const std::complex<double> * m_di, * m_gg, * m_rp;
    std::complex<double> * m_r, * m_x0;

    std::size_t m_m, m_m_curr;
    std::complex<double> * m_t, * m_w, ** m_VT, ** m_H, * m_d;
    preconditioners::preconditioner_interface<std::complex<double>, std::size_t> * m_precond;

    int m_num_threads;
    std::complex<double> * m_mv_tmp;
    std::size_t * m_mv_ind;
};

}}}}} // namespace core::solvers::CSRC::symmetric::complex

#endif // SOLVERS_CSRC_SYMMETRIC_COMPLEX_GMRES_OPENMP_H_INCLUDED
