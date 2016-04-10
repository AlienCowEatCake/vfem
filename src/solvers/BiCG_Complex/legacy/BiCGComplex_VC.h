#if !defined BICGCOMPLEX_VC_H_INCLUDED
#define BICGCOMPLEX_VC_H_INCLUDED

#include "../../../common/common.h"
#include "../../solver_interface.h"

#include <cstdlib>
#include <complex>
using namespace std;

class BiCGComplex_VC : public solver_interface<complex<double>, size_t>
{
public:
    void init(const size_t * gi_s, const size_t * gj_s, const complex<double> * di_s,
              const complex<double> * gg_s, const size_t n_s);
    void solve(complex<double> * solution, const complex<double> * rp, double gamma, size_t max_iter);

    BiCGComplex_VC();
    ~BiCGComplex_VC();
private:
    void mul_matrix(const complex<double> * f, complex<double> * x) const;
    complex<double> dot_prod(const complex<double> * a, const complex<double> * b) const;
    complex<double> dot_prod_nocj(const complex<double> * a, const complex<double> * b) const;
    double dot_prod_self(const complex<double> * a) const;

    size_t n;
    const size_t * gi, * gj;
    const complex<double> * di, * gg;

    complex<double> * r, * z, * s, * x0, * p, * t, * t1;
};

#endif // BICGCOMPLEX_VC_H_INCLUDED
