#if !defined BICG_COMPLEX_SMOOTH_H_INCLUDED
#define BICG_COMPLEX_SMOOTH_H_INCLUDED

#include "../../../common/common.h"
#include "../../solver_interface.h"

#include <cstdlib>
#include <complex>
using namespace std;

class BiCG_Complex_Smooth : public solver_interface<complex<double>, size_t>
{
public:
    void init(const size_t * gi_s, const size_t * gj_s, const complex<double> * di_s,
              const complex<double> * gg_s, const size_t n_s);
    void solve(complex<double> * solution, const complex<double> * rp, double gamma, size_t max_iter);

    BiCG_Complex_Smooth();
    ~BiCG_Complex_Smooth();
private:
    void mul_matrix(const complex<double> * f, complex<double> * x) const;
    complex<double> dot_prod(const complex<double> * a, const complex<double> * b) const;
    complex<double> dot_prod_nocj(const complex<double> * a, const complex<double> * b) const;
    double dot_prod_self(const complex<double> * a) const;
    double dot_prod_real(const complex<double> * a, const complex<double> * b) const;

    size_t n;
    const size_t * gi, * gj;
    const complex<double> * di, * gg;

    complex<double> * r, * z, * s, * x0, * p, * t, * t1, * xs, * rs;
};

#endif // BICG_COMPLEX_SMOOTH_H_INCLUDED
