#ifndef GAUSS_H_INCLUDED
#define GAUSS_H_INCLUDED

#include <cmath>
#include <cstdlib>
#include <complex>

using namespace std;

class Gauss
{
public:
    // Инициализация несимметричная
    void init(size_t * gi_s, size_t * gj_s, complex<double> * di_s,
              complex<double> * gl_s, complex<double> * gu_s, size_t n_s);
    // Инициализация симметричная
    void init(size_t * gi_s, size_t * gj_s, complex<double> * di_s,
              complex<double> * gg_s, size_t n_s);
    // Получение решения
    void solve(complex<double> * solution, complex<double> * rp_s, double eps, size_t max_iter);

    Gauss();
    ~Gauss();
private:
    size_t n_gauss;
    complex<double> ** A;

    void convert_matrix();
    void delete_matrix();

    size_t n;
    size_t * gi, * gj;
    complex<double> * di, * gl, * gu, * rp, * r, * x0, * z, * p, * s, * xs, * rs;
    complex<double> * L_di, * L_gg;
};


#endif // GAUSS_H_INCLUDED
