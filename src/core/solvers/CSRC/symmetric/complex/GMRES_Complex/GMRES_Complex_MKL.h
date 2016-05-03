#if !defined(SOLVERS_CSRC_SYMMETRIC_COMPLEX_GMRES_MKL_H_INCLUDED)
#define SOLVERS_CSRC_SYMMETRIC_COMPLEX_GMRES_MKL_H_INCLUDED

#include <cstdlib>
#include <complex>
#include "../../symmetric_solver_interface.h"
#include "../../../../../wrappers/mkl_wrapper.h"

namespace core { namespace solvers { namespace CSRC { namespace symmetric { namespace complex {

/**
 * @brief Комплексный GMRES, использующий библиотеку MKL
 */
class GMRES_Complex_MKL : public symmetric_solver_interface<std::complex<double>, std::size_t>
{
public:
    void init(const std::size_t * gi, const std::size_t * gj, const std::complex<double> * di,
              const std::complex<double> * gg, std::size_t n,
              preconditioners::preconditioner_interface<std::complex<double>, std::size_t> * precond);
    void solve(std::complex<double> * solution, const std::complex<double> * rp,
               double eps, std::size_t max_iter);

    GMRES_Complex_MKL();
    ~GMRES_Complex_MKL();

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
    MKL_INT m_m_mkl;
    MKL_INT * m_ia, * m_ja;
    MKL_Complex16 * m_aa;
};

}}}}} // namespace core::solvers::CSRC::symmetric::complex

#endif // SOLVERS_CSRC_SYMMETRIC_COMPLEX_GMRES_MKL_H_INCLUDED
