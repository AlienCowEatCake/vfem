#if !defined COCG_LLT_SMOOTH_H_INCLUDED
#define COCG_LLT_SMOOTH_H_INCLUDED

#include <cstdlib>
#include <complex>
using namespace std;

class COCG_LLT_Smooth
{
public:
    void init(size_t * gi_s, size_t * gj_s, complex<double> * di_s,
              complex<double> * gg_s, size_t n_s);
    void solve(complex<double> * solution, complex<double> * rp_s, double eps);

    COCG_LLT_Smooth();
    ~COCG_LLT_Smooth();
private:
    void make_LLT_decomposition();
    void mul_matrix(const complex<double> * f, complex<double> * x) const;
    void solve_L(const complex<double> * f, complex<double> * x) const;
    void solve_LT(complex<double> * f, complex<double> * x) const;
    void solve_LLT(const complex<double> * f, complex<double> * x) const;
    complex<double> dot_prod_nocj(const complex<double> * a, const complex<double> * b) const;
    double dot_prod_self(const complex<double> * a) const;
    double dot_prod_real(const complex<double> * a, const complex<double> * b) const;
    bool is_fpu_error(double x) const;

    size_t n;
    size_t * gi, * gj;
    complex<double> * di, * gg, * rp, * r, * x0, * z, * p, * s, * xs, * rs;
    complex<double> * L_di, * L_gg;
};

#endif // COCG_LLT_SMOOTH_H_INCLUDED
