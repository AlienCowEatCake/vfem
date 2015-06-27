#include "slae.h"

// ============================================================================

SLAE::SLAE()
{
    gg = di = rp = x = NULL;
    ig = jg = NULL;
    n = 0;
}

SLAE::~SLAE()
{
    dealloc_all();
}

void SLAE::solve(double eps)
{
    cout << "Solving SLAE ..." << endl;
    solver.init(ig, jg, di, gg, n);
    solver.solve(x, rp, eps);
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

void SLAE::dump(const string & filename) const
{
    ofstream slae_file;
    slae_file.open(filename.c_str(), ios::out);

    if(!slae_file.good())
    {
        cerr << "Error in " << __FILE__ << ":" << __LINE__
             << " while writing file " << filename << endl;
        throw IO_FILE_ERROR;
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
}

void SLAE::restore(const string & filename)
{
    ifstream slae_file;
    slae_file.open(filename.c_str(), ios::in);
    if(!slae_file.good())
    {
        cerr << "Error in " << __FILE__ << ":" << __LINE__
             << " while reading file " << filename << endl;
        throw IO_FILE_ERROR;
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
        cerr << "Error in " << __FILE__ << ":" << __LINE__
             << " while reading file " << filename << endl;
        throw IO_FILE_ERROR;
    }

    slae_file.close();
}

void SLAE::dump_x(const string & filename) const
{
    ofstream slae_file;
    slae_file.open(filename.c_str(), ios::out);

    if(!slae_file.good())
    {
        cerr << "Error in " << __FILE__ << ":" << __LINE__
             << " while writing file " << filename << endl;
        throw IO_FILE_ERROR;
    }

    slae_file.precision(17);
    slae_file.setf(ios::scientific);

    for(size_t i = 0; i < n; i++)
        slae_file << x[i] << '\n';

    slae_file << '\n';
    slae_file.flush();
    slae_file.close();
}

void SLAE::restore_x(const string & filename)
{
    ifstream slae_file;
    slae_file.open(filename.c_str(), ios::in);
    if(!slae_file.good())
    {
        cerr << "Error in " << __FILE__ << ":" << __LINE__
             << " while reading file " << filename << endl;
        throw IO_FILE_ERROR;
    }

    for(size_t i = 0; i < n; i++)
        slae_file >> x[i];

    if(!slae_file.good())
    {
        cerr << "Error in " << __FILE__ << ":" << __LINE__
             << " while reading file " << filename << endl;
        throw IO_FILE_ERROR;
    }

    slae_file.close();
}

// ============================================================================

SLAE_V_cycle::SLAE_V_cycle()
{
    ker_gg = ker_di = NULL;
    ker_ig = ker_jg = NULL;
    ker_n = proj_n = 0;
    proj_gg = NULL;
    proj_ig = proj_jg = NULL;
}

SLAE_V_cycle::~SLAE_V_cycle()
{
    if(ker_gg) delete [] ker_gg;
    if(ker_di) delete [] ker_di;
    if(ker_ig) delete [] ker_ig;
    if(ker_jg) delete [] ker_jg;
    if(proj_gg) delete [] proj_gg;
    if(proj_ig) delete [] proj_ig;
    if(proj_jg) delete [] proj_jg;
}

void SLAE_V_cycle::solve(double eps, set<size_t> * edges_first, set<size_t> * ker_edges_first)
{
    cout << "Solving V-Cycle SLAE ..." << endl;
    cout << "Main dim: " << n << "\nKernel dim: " << ker_n << endl;

    solver_V.init_lvl1(ig, jg, di, gg, rp, n);
    solver_V.init_lvl2(ker_ig, ker_jg, ker_di, ker_gg, ker_n);
    solver_V.get_main_bounds(edges_first);
    solver_V.get_grad_bounds(ker_edges_first);
    solver_V.init_operators(proj_ig, proj_jg, proj_gg);

    solver_V.solve(x, eps);
}

void SLAE_V_cycle::ker_alloc_all(size_t n_size, size_t gg_size)
{
    ker_n = n_size;
    ker_ig = new size_t [ker_n + 1];
    ker_jg = new size_t [gg_size];
    ker_di = new complex<double> [ker_n];
    ker_gg = new complex<double> [gg_size];

    for(size_t i = 0; i < ker_n; i++)
        ker_di[i] = complex<double>(0.0, 0.0);
    for(size_t i = 0; i < gg_size; i++)
        ker_gg[i] = complex<double>(0.0, 0.0);
}

void SLAE_V_cycle::ker_add(size_t i, size_t j, complex<double> elem)
{
    assert(i < ker_n);
    assert(j < ker_n);
    if(j > i)
        swap(i, j);
    size_t ind = 0;
    bool flag = false;
    for(size_t k = ker_ig[i]; k < ker_ig[i + 1] && !flag; k++)
    {
        if(ker_jg[k] == j)
        {
            ind = k;
            flag = true;
        }
    }
    assert(flag != false);
    ker_gg[ind] += elem;
}

void SLAE_V_cycle::proj_alloc_all(size_t n_size, size_t gg_size)
{
    proj_n = n_size;
    proj_ig = new size_t [proj_n + 1];
    proj_jg = new size_t [gg_size];
    proj_gg = new double [gg_size];

    for(size_t i = 0; i < gg_size; i++)
        proj_gg[i] = 0.0;
}

void SLAE_V_cycle::proj_set(size_t i, size_t j, double elem)
{
    assert(i < ker_n);
    assert(j < n);
    size_t ind = 0;
    bool flag = false;
    for(size_t k = proj_ig[i]; k < proj_ig[i + 1] && !flag; k++)
    {
        if(proj_jg[k] == j)
        {
            ind = k;
            flag = true;
        }
    }
    assert(flag != false);
    proj_gg[ind] = elem;
}
