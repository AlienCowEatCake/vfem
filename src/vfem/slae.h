#if !defined SLAE_H_INCLUDED
#define SLAE_H_INCLUDED

#include "../common/common.h"
//#include "../solvers/BiCGComplex_VC.h"
//#include "../solvers/BiCGStabComplex_VC.h"
//#include "../solvers/CGMComplex_VC.h"
//#include "../solvers/CGMComplex_LLT.h"
#include "../solvers/COCG_LLT_Smooth.h"
#include "../solvers/COCG_LLT_Smooth_MKL.h"

// Класс СЛАУ
class SLAE
{
public:
    SLAE();
    ~SLAE();
    void solve(double eps);
    inline void inline_init()
    {
        solver.init(ig, jg, di, gg, n);
    }
    inline void inline_solve(complex<double> * solution, complex<double> * rp_s, double eps, size_t max_iter)
    {
        solver.solve(solution, rp_s, eps, max_iter);
    }
    void alloc_all(size_t n_size, size_t gg_size);
    void dealloc_all();
    void add(size_t i, size_t j, const complex<double> & elem);
    complex<double> * gg, * di, * rp, * x;
    size_t * ig, * jg;
    size_t n;
    bool dump(const string & filename) const;
    bool restore(const string & filename);
    bool dump_x(const string & filename) const;
    bool restore_x(const string & filename);
private:
    //BiCGComplex_VC solver;
    //BiCGStabComplex_VC solver;
    //CGMComplex_VC solver;
    //CGMComplex_LLT solver;
#if !defined USE_MKL
    COCG_LLT_Smooth solver;
#else
    COCG_LLT_Smooth_MKL solver;
#endif
};

#endif // SLAE_H_INCLUDED
