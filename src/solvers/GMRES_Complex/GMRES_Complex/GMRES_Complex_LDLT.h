#if !defined(GMRES_COMPLEX_LDLT_H_INCLUDED)
#define GMRES_COMPLEX_LDLT_H_INCLUDED

#include <cstdlib>
#include <complex>
#include "../../solver_interface.h"

using namespace std;

class GMRES_Complex_LDLT : public solver_interface<complex<double>, size_t>
{
public:
    void init(const size_t * gi_s, const size_t * gj_s, const complex<double> * di_s,
              const complex<double> * gg_s, size_t n_s);
    void solve(complex<double> * solution, const complex<double> * rp_s,
               double eps, size_t max_iter);

    GMRES_Complex_LDLT();
    ~GMRES_Complex_LDLT();
private:
    void make_LDLT_decomposition();
    void mul_matrix(const complex<double> * f, complex<double> * x) const;
    void solve_LD(const complex<double> * f, complex<double> * x) const;
    void solve_LT(complex<double> * f, complex<double> * x) const;
    void mul_LT(const complex<double> * f, complex<double> * x) const;
    void solve_LTALD(const complex<double> * f, complex<double> * x, complex<double> * tmp) const;
    complex<double> dot_prod(const complex<double> * a, const complex<double> * b) const;
    double dot_prod_self(const complex<double> * a) const;
    bool is_fpu_error(double x) const;
    void copy_vec(const complex<double> * f, complex<double> * x) const;

    size_t n;
    const size_t * gi, * gj;
    const complex<double> * di, * gg, * rp;
    complex<double> * r, * x0;
    complex<double> * L_di, * L_gg;

    size_t m, m_curr;
    complex<double> * t, * w, ** VT, ** H, * d;
};

#endif // GMRES_COMPLEX_LDLT_H_INCLUDED
