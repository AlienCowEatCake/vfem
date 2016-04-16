#include "slae.h"

// Решатели COCG
#include "../solvers/COCG/COCG/COCG.h"
#include "../solvers/COCG/COCG/COCG_Di.h"
#include "../solvers/COCG/COCG/COCG_GS.h"
#include "../solvers/COCG/COCG/COCG_LLT.h"
#include "../solvers/COCG/COCG/COCG_LDLT.h"
#include "../solvers/COCG/COCG/COCG_Smooth.h"
#include "../solvers/COCG/COCG/COCG_Di_Smooth.h"
#include "../solvers/COCG/COCG/COCG_GS_Smooth.h"
#include "../solvers/COCG/COCG/COCG_LLT_Smooth.h"
#include "../solvers/COCG/COCG/COCG_LDLT_Smooth.h"

// Решатели COCG с OpenMP
#include "../solvers/COCG/COCG_OpenMP/COCG_Smooth_OpenMP.h"
#include "../solvers/COCG/COCG_OpenMP/COCG_Di_Smooth_OpenMP.h"

// Решатели COCG с MKL
#include "../solvers/COCG/COCG_MKL/COCG_Smooth_MKL.h"
#include "../solvers/COCG/COCG_MKL/COCG_Di_Smooth_MKL.h"
#include "../solvers/COCG/COCG_MKL/COCG_LLT_Smooth_MKL.h"

// Решатели GMRES
#include "../solvers/GMRES_Complex/GMRES_Complex/GMRES_Complex.h"
#include "../solvers/GMRES_Complex/GMRES_Complex/GMRES_Complex_Di.h"
#include "../solvers/GMRES_Complex/GMRES_Complex/GMRES_Complex_LLT.h"
#include "../solvers/GMRES_Complex/GMRES_Complex/GMRES_Complex_LDLT.h"

// Решатели GMRES с OpenMP
#include "../solvers/GMRES_Complex/GMRES_Complex_OpenMP/GMRES_Complex_OpenMP.h"
#include "../solvers/GMRES_Complex/GMRES_Complex_OpenMP/GMRES_Complex_Di_OpenMP.h"

// Решатели GMRES с MKL
#include "../solvers/GMRES_Complex/GMRES_Complex_MKL/GMRES_Complex_MKL.h"
#include "../solvers/GMRES_Complex/GMRES_Complex_MKL/GMRES_Complex_Di_MKL.h"
#include "../solvers/GMRES_Complex/GMRES_Complex_MKL/GMRES_Complex_LLT_MKL.h"

// Решатели BiCG
#include "../solvers/BiCG_Complex/BiCG_Complex/BiCG_Complex_Smooth.h"

// Решатели BiCGStab
#include "../solvers/BiCGStab_Complex/BiCGStab_Complex/BiCGStab_Complex_Smooth.h"

// Решатели COCR
#include "../solvers/COCR/COCR/COCR.h"
#include "../solvers/COCR/COCR/COCR_Di.h"
#include "../solvers/COCR/COCR/COCR_Di_Smooth.h"
#include "../solvers/COCR/COCR/COCR_LLT.h"
#include "../solvers/COCR/COCR/COCR_LLT_Smooth.h"
#include "../solvers/COCR/COCR/COCR_LDLT.h"
#include "../solvers/COCR/COCR/COCR_LDLT_Smooth.h"

// Legacy решатели
#include "../solvers/BiCG_Complex/legacy/BiCGComplex_VC.h"
#include "../solvers/BiCGStab_Complex/legacy/BiCGStabComplex_VC.h"
#include "../solvers/COCG/legacy/CGMComplex_VC.h"
#include "../solvers/COCG/legacy/CGMComplex_LLT.h"

SLAE::SLAE()
{
    gg = di = rp = x = NULL;
    ig = jg = NULL;
    n = 0;
    solver = NULL;
}

SLAE::~SLAE()
{
    dealloc_all();
    delete solver;
}

void SLAE::init(const string & name)
{
    delete solver;
    // Решатели COCG
    if(name == "COCG")
        solver = new COCG;
    else if(name == "COCG_Di")
        solver = new COCG_Di;
    else if(name == "COCG_GS")
        solver = new COCG_GS;
    else if(name == "COCG_LLT")
        solver = new COCG_LLT;
    else if(name == "COCG_LDLT")
        solver = new COCG_LDLT;
    else if(name == "COCG_Smooth")
        solver = new COCG_Smooth;
    else if(name == "COCG_Di_Smooth")
        solver = new COCG_Di_Smooth;
    else if(name == "COCG_GS_Smooth")
        solver = new COCG_GS_Smooth;
    else if(name == "COCG_LLT_Smooth")
        solver = new COCG_LLT_Smooth;
    else if(name == "COCG_LDLT_Smooth")
        solver = new COCG_LDLT_Smooth;
    // Решатели COCG с OpenMP
    else if(name == "COCG_Smooth_OpenMP")
        solver = new COCG_Smooth_OpenMP;
    else if(name == "COCG_Di_Smooth_OpenMP")
        solver = new COCG_Di_Smooth_OpenMP;
    // Решатели COCG с MKL
    else if(name == "COCG_Smooth_MKL")
        solver = new COCG_Smooth_MKL;
    else if(name == "COCG_Di_Smooth_MKL")
        solver = new COCG_Di_Smooth_MKL;
    else if(name == "COCG_LLT_Smooth_MKL")
        solver = new COCG_LLT_Smooth_MKL;
    // Решатели GMRES
    else if(name == "GMRES_Complex")
        solver = new GMRES_Complex;
    else if(name == "GMRES_Complex_Di")
        solver = new GMRES_Complex_Di;
    else if(name == "GMRES_Complex_LLT")
        solver = new GMRES_Complex_LLT;
    else if(name == "GMRES_Complex_LDLT")
        solver = new GMRES_Complex_LDLT;
    // Решатели GMRES с OpenMP
    else if(name == "GMRES_Complex_OpenMP")
        solver = new GMRES_Complex_OpenMP;
    else if(name == "GMRES_Complex_Di_OpenMP")
        solver = new GMRES_Complex_Di_OpenMP;
    // Решатели GMRES с MKL
    else if(name == "GMRES_Complex_MKL")
        solver = new GMRES_Complex_MKL;
    else if(name == "GMRES_Complex_Di_MKL")
        solver = new GMRES_Complex_Di_MKL;
    else if(name == "GMRES_Complex_LLT_MKL")
        solver = new GMRES_Complex_LLT_MKL;
    // Решатели BiCG
    else if(name == "BiCG_Complex_Smooth")
        solver = new BiCG_Complex_Smooth;
    // Решатели BiCGStab
    else if(name == "BiCGStab_Complex_Smooth")
        solver = new BiCGStab_Complex_Smooth;
    // Решатели COCR
    else if(name == "COCR")
        solver = new COCR;
    else if(name == "COCR_Di")
        solver = new COCR_Di;
    else if(name == "COCR_Di_Smooth")
        solver = new COCR_Di_Smooth;
    else if(name == "COCR_LLT")
        solver = new COCR_LLT;
    else if(name == "COCR_LLT_Smooth")
        solver = new COCR_LLT_Smooth;
    else if(name == "COCR_LDLT")
        solver = new COCR_LDLT;
    else if(name == "COCR_LDLT_Smooth")
        solver = new COCR_LDLT_Smooth;
    // Legacy решатели
    else if(name == "BiCGComplex_VC")
        solver = new BiCGComplex_VC;
    else if(name == "BiCGStabComplex_VC")
        solver = new BiCGStabComplex_VC;
    else if(name == "CGMComplex_VC")
        solver = new CGMComplex_VC;
    else if(name == "CGMComplex_LLT")
        solver = new CGMComplex_LLT;
    // Решатель по умолчанию
    else
    {
        cout << "Warning in " << __FILE__ << ":" << __LINE__
             << " - unknown solver \"" << name << "\"" << endl;
#if !defined(USE_MKL)
        solver = new COCG_LLT_Smooth;
#else
        solver = new COCG_LLT_Smooth_MKL;
#endif
    }
    solver->init(ig, jg, di, gg, n);
}

void SLAE::solve(const string & name, double eps, size_t max_iter)
{
    cout << "Solving SLAE ..." << endl;
    init(name);
    step_solve(x, rp, eps, max_iter);
}

void SLAE::alloc_all(size_t n_size, size_t gg_size)
{
    n = n_size;
    ig = new size_t [n + 1];
    jg = new size_t [gg_size];
    di = new complex<double> [n];
    gg = new complex<double> [gg_size];
    rp = new complex<double> [n];
    x  = new complex<double> [n];

    for(size_t i = 0; i < n; i++)
    {
        di[i] = complex<double>(0.0, 0.0);
        rp[i] = complex<double>(0.0, 0.0);
        x[i]  = complex<double>(0.0, 0.0);
    }
    for(size_t i = 0; i < gg_size; i++)
    {
        gg[i] = complex<double>(0.0, 0.0);
    }
}

void SLAE::dealloc_all()
{
    delete [] gg;
    delete [] di;
    delete [] rp;
    delete [] x;
    delete [] ig;
    delete [] jg;
    gg = di = rp = x = NULL;
    ig = jg = NULL;
    n = 0;
}

void SLAE::add(size_t i, size_t j, const complex<double> & elem)
{
    assert(i < n);
    assert(j < n);
    if(j > i)
        swap(i, j);
    size_t ind = 0;
    bool flag = false;
    for(size_t k = ig[i]; k < ig[i + 1] && !flag; k++)
    {
        if(jg[k] == j)
        {
            ind = k;
            flag = true;
        }
    }
    assert(flag != false);
    gg[ind] += elem;
}

bool SLAE::dump(const string & filename) const
{
    ofstream slae_file;
    slae_file.open(filename.c_str(), ios::out);

    if(!slae_file.good())
    {
        cout << "Error in " << __FILE__ << ":" << __LINE__
             << " while writing file " << filename << endl;
        return false;
    }

    slae_file.precision(17);
    slae_file.setf(ios::scientific);

    size_t gg_size = ig[n];
    slae_file << n << ' ' << gg_size << '\n';

    for(size_t i = 0; i <= n; i++)
        slae_file << ig[i] << '\n';

    for(size_t i = 0; i < gg_size; i++)
        slae_file << jg[i] << '\n';

    for(size_t i = 0; i < n; i++)
        slae_file << di[i] << ' ' << rp[i] << ' ' << x[i] << '\n';

    for(size_t i = 0; i < gg_size; i++)
        slae_file << gg[i] << '\n';

    slae_file << '\n';
    slae_file.flush();
    slae_file.close();
    return true;
}

bool SLAE::restore(const string & filename)
{
    ifstream slae_file;
    slae_file.open(filename.c_str(), ios::in);
    if(!slae_file.good())
    {
        cout << "Error in " << __FILE__ << ":" << __LINE__
             << " while reading file " << filename << endl;
        return false;
    }

    size_t gg_size;
    dealloc_all();
    slae_file >> n >> gg_size;
    alloc_all(n, gg_size);

    for(size_t i = 0; i <= n; i++)
        slae_file >> ig[i];

    for(size_t i = 0; i < gg_size; i++)
        slae_file >> jg[i];

    for(size_t i = 0; i < n; i++)
        slae_file >> di[i] >> rp[i] >> x[i];

    for(size_t i = 0; i < gg_size; i++)
        slae_file >> gg[i];

    if(!slae_file.good())
    {
        cout << "Error in " << __FILE__ << ":" << __LINE__
             << " while reading file " << filename << endl;
        return false;
    }

    slae_file.close();
    return true;
}

bool SLAE::dump_x(const string & filename) const
{
    ofstream slae_file;
    slae_file.open(filename.c_str(), ios::out);

    if(!slae_file.good())
    {
        cout << "Error in " << __FILE__ << ":" << __LINE__
             << " while writing file " << filename << endl;
        return false;
    }

    slae_file.precision(17);
    slae_file.setf(ios::scientific);

    for(size_t i = 0; i < n; i++)
        slae_file << x[i] << '\n';

    slae_file << '\n';
    slae_file.flush();
    slae_file.close();
    return true;
}

bool SLAE::restore_x(const string & filename)
{
    ifstream slae_file;
    slae_file.open(filename.c_str(), ios::in);
    if(!slae_file.good())
    {
        cout << "Error in " << __FILE__ << ":" << __LINE__
             << " while reading file " << filename << endl;
        return false;
    }

    for(size_t i = 0; i < n; i++)
        slae_file >> x[i];

    if(!slae_file.good())
    {
        cout << "Error in " << __FILE__ << ":" << __LINE__
             << " while reading file " << filename << endl;
        return false;
    }

    slae_file.close();
    return true;
}
