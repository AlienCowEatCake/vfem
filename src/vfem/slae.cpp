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

void SLAE::solve(double eps, size_t max_iter)
{
    cout << "Solving SLAE ..." << endl;
    solver.init(ig, jg, di, gg, n);
    solver.solve(x, rp, eps, max_iter);
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

    if(i == j)
    {
        di[i] += elem;
        return;
    }

    if(j > i)
    {
        return;
        swap(i, j);
    }
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

// ============================================================================

SLAE_ns::SLAE_ns()
{
    gl = gu = di = rp = x = NULL;
    ig = jg = NULL;
    n = 0;
}

SLAE_ns::~SLAE_ns()
{
    dealloc_all();
}

void SLAE_ns::solve(double eps, size_t max_iter)
{
    cout << "Solving SLAE ..." << endl;
    solver.init(ig, jg, di, gl, gu, n);
    solver.solve(x, rp, eps, max_iter);
}

void SLAE_ns::alloc_all(size_t n_size, size_t gg_size)
{
    n = n_size;
    ig = new size_t [n + 1];
    jg = new size_t [gg_size];
    di = new complex<double> [n];
    gl = new complex<double> [gg_size];
    gu = new complex<double> [gg_size];
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
        gl[i] = complex<double>(0.0, 0.0);
        gu[i] = complex<double>(0.0, 0.0);
    }
}

void SLAE_ns::dealloc_all()
{
    delete [] gl;
    delete [] gu;
    delete [] di;
    delete [] rp;
    delete [] x;
    delete [] ig;
    delete [] jg;
    gl = gu = di = rp = x = NULL;
    ig = jg = NULL;
    n = 0;
}

void SLAE_ns::add(size_t i, size_t j, const complex<double> & elem)
{
    assert(i < n);
    assert(j < n);

    if(i == j)
    {
        di[i] += elem;
        return;
    }

    size_t i_p = i, j_p = j;
    if(j_p > i_p)
        swap(i_p, j_p);
    size_t ind = 0;
    bool flag = false;
    for(size_t k = ig[i_p]; k < ig[i_p + 1] && !flag; k++)
    {
        if(jg[k] == j_p)
        {
            ind = k;
            flag = true;
        }
    }

    if(j < i)
        gl[ind] += elem;
    else
        gu[ind] += elem;
}

bool SLAE_ns::dump(const string & filename) const
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
        slae_file << gl[i] << ' ' << gu[i] << '\n';

    slae_file << '\n';
    slae_file.flush();
    slae_file.close();
    return true;
}

bool SLAE_ns::restore(const string & filename)
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
        slae_file >> gl[i] >> gu[i];

    if(!slae_file.good())
    {
        cout << "Error in " << __FILE__ << ":" << __LINE__
             << " while reading file " << filename << endl;
        return false;
    }

    slae_file.close();
    return true;
}

bool SLAE_ns::dump_x(const string & filename) const
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

bool SLAE_ns::restore_x(const string & filename)
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

// ============================================================================
