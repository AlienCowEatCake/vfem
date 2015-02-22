#if !defined BICGCOMPLEX_VC_H_INCLUDED
#define BICGCOMPLEX_VC_H_INCLUDED

#include "../common/common.h"

#include <cstdlib>
#include <complex>
using namespace std;

class BiCGComplex_VC
{
public:
    void init(size_t * gi_s, size_t * gj_s, complex<double> * di_s,
              complex<double> * gg_s, size_t n_s);
    void solve(complex<double> * solution, complex<double> * rp, double gamma);

    BiCGComplex_VC();
    ~BiCGComplex_VC();
private:
    void mul_matrix(const complex<double> * f, complex<double> * x) const;
    complex<double> dot_prod(const complex<double> * a, const complex<double> * b) const;
    complex<double> dot_prod_nocj(const complex<double> * a, const complex<double> * b) const;

    size_t n;
    size_t * gi, * gj;
    complex<double> * di, * gg;

    complex<double> * r, * z, * s, * x0, * p, * t, * t1;
};

#endif // BICGCOMPLEX_VC_H_INCLUDED
