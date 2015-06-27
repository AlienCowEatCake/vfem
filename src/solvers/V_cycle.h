#if !defined V_CYCLE_H_INCLUDED
#define V_CYCLE_H_INCLUDED

#include "CGMComplex_VC.h"
#include "BiCGComplex_VC.h"
#include "BiCGStabComplex_VC.h"
#include "COCG_LLT_Smooth.h"

#include <cstdlib>
#include <complex>
#include <set>
using namespace std;

/*
== Двухуровневый решатель V-цикл для решения гармонического векторного уравнения Гельмгольца ==

	Основное решение идёт на полном базисе, предобусловлевание на градиентном базисе.
	Градиентный базис состоит из градиентов скалярных (узловых) базисных функций.
	Приставка lvl1 означает первый уровень (мелкая сетка) - роторные функции
	Приставка lvl2 означает второй уровень (грубая сетка) - градиентные функции
	В качестве базового решатель используется BiCGStub.
*/

class V_cycle
{
public:
    void init_lvl1(size_t * gi_s, size_t * gj_s, complex<double> * di_s, complex<double> * gg_s, complex<double> * rp_s, size_t n_h_s);
    void init_lvl2(size_t * gi_s, size_t * gj_s, complex<double> * di_s, complex<double> * gg_s, size_t n_2h_s);
    void get_grad_bounds(set<size_t> * grad_bounds_s);
    void get_main_bounds(set<size_t> * main_bounds_s);

    void init_operators(size_t * gi_R_s, size_t * gj_R_s, double * gg_R_s);
    void solve(complex<double> * solution, double eps);

private:
    void mul_matrix(const complex<double> * f, complex<double> * x) const;
    complex<double> dot_prod(const complex<double> * a, const complex<double> * b) const;
    void calc_residual(const complex<double> * x0, complex<double> * p) const;

    size_t n_lvl1, n_lvl2;

    size_t * gi_lvl1, * gj_lvl1;
    complex<double> * di_lvl1, * gg_lvl1, * rp;

    size_t * gi_lvl2, * gj_lvl2;
    complex<double> * di_lvl2, * gg_lvl2;

    size_t * gi_R, * gj_R;
    double * gg_R;

    set<size_t> * grad_bounds;
    set<size_t> * main_bounds;

    //CGMComplex_VC bcgm_lvl1;
    //CGMComplex_VC bcgm_lvl2;
    COCG_LLT_Smooth bcgm_lvl1;
    COCG_LLT_Smooth bcgm_lvl2;
    //BiCGStabComplex_VC bcgm_lvl1;
    //BiCGStabComplex_VC bcgm_lvl2;
};

#endif
