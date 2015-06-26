#if !defined SLAE_H_INCLUDED
#define SLAE_H_INCLUDED

#include "../common/common.h"
//#include "../solvers/BiCGComplex_VC.h"
//#include "../solvers/BiCGStabComplex_VC.h"
//#include "../solvers/CGMComplex_VC.h"
//#include "../solvers/CGMComplex_LLT.h"
#include "../solvers/COCG_LLT_Smooth.h"

#include "../solvers/V_cycle.h"

// Класс СЛАУ
class SLAE
{
public:
    SLAE();
    ~SLAE();
    void solve(double eps);
    void alloc_all(size_t n_size, size_t gg_size);
    void dealloc_all();
    void add(size_t i, size_t j, const complex<double> & elem);
    complex<double> * gg, * di, * rp, * x;
    size_t * ig, * jg;
    size_t n;
    void dump(const string & filename) const;
    void restore(const string & filename);
    void dump_x(const string & filename) const;
    void restore_x(const string & filename);
private:
    //BiCGComplex_VC solver;
    //BiCGStabComplex_VC solver;
    //CGMComplex_VC solver;
    //CGMComplex_LLT solver;
    COCG_LLT_Smooth solver;
};

// Класс СЛАУ для V-цикла
class SLAE_V_cycle : public SLAE
{
public:
    SLAE_V_cycle();
    ~SLAE_V_cycle();
    void solve(double eps, set<size_t> * edges_first, set<size_t> * ker_edges_first);
    // Ядро
    void ker_alloc_all(size_t n_size, size_t gg_size);
    void ker_add(size_t i, size_t j, complex<double> elem);
    complex<double> * ker_gg, * ker_di;
    size_t * ker_ig, * ker_jg;
    size_t ker_n;
    // Проектор
    void proj_alloc_all(size_t n_size, size_t gg_size);
    void proj_set(size_t i, size_t j, double elem);
    double * proj_gg;
    size_t * proj_ig, * proj_jg;
    size_t proj_n;
private:
    V_cycle solver_V;
};

#endif // SLAE_H_INCLUDED
