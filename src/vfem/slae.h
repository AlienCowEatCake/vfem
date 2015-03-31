#if !defined SLAE_H_INCLUDED
#define SLAE_H_INCLUDED

#include "../common/common.h"
//#include "../solvers/BiCGComplex_VC.h"
//#include "../solvers/BiCGStabComplex_VC.h"
//#include "../solvers/CGMComplex_VC.h"
//#include "../solvers/CGMComplex_LLT.h"
#include "../solvers/COCG_LLT_Smooth.h"

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

#endif // SLAE_H_INCLUDED
