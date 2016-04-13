#if !defined(COCG_LLT_H_INCLUDED)
#define COCG_LLT_H_INCLUDED

#include <cstdlib>
#include <complex>
#include "../../solver_interface.h"

using namespace std;

class COCG_LLT : public solver_interface<complex<double>, size_t>
{
public:
    void init(const size_t * gi_s, const size_t * gj_s, const complex<double> * di_s,
              const complex<double> * gg_s, size_t n_s);
    void solve(complex<double> * solution, const complex<double> * rp_s, double eps, size_t max_iter);

    COCG_LLT();
    ~COCG_LLT();
private:
    void make_LLT_decomposition();
    void mul_matrix(const complex<double> * f, complex<double> * x) const;
    void solve_L(const complex<double> * f, complex<double> * x) const;
    void solve_LT(complex<double> * f, complex<double> * x) const;
    void solve_LLT(const complex<double> * f, complex<double> * x) const;
    complex<double> dot_prod_nocj(const complex<double> * a, const complex<double> * b) const;
    double dot_prod_self(const complex<double> * a) const;
    bool is_fpu_error(double x) const;

    size_t n;
    const size_t * gi, * gj;
    const complex<double> * di, * gg, * rp;
    complex<double> * r, * x0, * z, * p, * s;
    complex<double> * L_di, * L_gg;
};

#endif // COCG_LLT_H_INCLUDED
