#include "problems.h"

#if defined ANALYTICAL_CUBE
#if !defined VFEM_USE_NONHOMOGENEOUS_FIRST || !defined VFEM_USE_ANALYTICAL || defined VFEM_USE_PML
#error "Please, reconfigure!"
#endif

double SLAE_SURF_EPSILON = 1e-15;
double SLAE_MAIN_EPSILON = 1e-15;

cvector3 func_true(const point & p)
{
    MAYBE_UNUSED(p);

    //return cvector3(1.0, 0.0, 0.0);
    //return cvector3(p.y + p.z, p.x + p.z, p.x + p.y);
    return cvector3(exp(p.y + p.z), exp(p.x + p.z), exp(p.x + p.y));
}

cvector3 func_rp(const point & p, const phys_area & phys)
{
    MAYBE_UNUSED(p);
    MAYBE_UNUSED(phys);

    complex<double> k2(- phys.omega * phys.omega * phys.epsilon, phys.omega * phys.sigma);
    //return k2 * cvector3(1.0, 0.0, 0.0);
    //return k2 * cvector3(- p.y - p.z, - p.x - p.z, - p.x - p.y);
    return cvector3(
               -2.0 * exp(p.y + p.z) / phys.mu + k2 * exp(p.y + p.z),
               -2.0 * exp(p.x + p.z) / phys.mu + k2 * exp(p.x + p.z),
               -2.0 * exp(p.x + p.y) / phys.mu + k2 * exp(p.x + p.y)
           );
}

cvector3 func_b1(const point & p, const triangle * tr)
{
    MAYBE_UNUSED(p);
    MAYBE_UNUSED(tr);

    cvector3 t = func_true(p);
    size_t phys_num = tr->phys->gmsh_num;
    switch(phys_num)
    {
    case 30:
        return cvector3(0.0, t.y, t.z);
    case 29:
        return cvector3(0.0, t.y, t.z);
    case 32:
        return cvector3(t.x, 0.0, t.z);
    case 33:
        return cvector3(t.x, 0.0, t.z);
    case 28:
        return cvector3(t.x, t.y, 0.0);
    case 31:
        return cvector3(t.x, t.y, 0.0);
    }
    cerr << "Unknown bound!" << endl;
    return cvector3(0.0, 0.0, 0.0);
}

//string mesh_filename = "data/cube/cube.msh";
//string phys_filename = "data/cube/phys_analyt.txt";
string mesh_filename = "data/approx/0.0250/cube.msh";
string phys_filename = "data/approx/phys_param.txt";
string tecplot_filename = "analytical_cube.plt";

void postprocessing(VFEM & v, char * timebuf)
{
    MAYBE_UNUSED(v);
    MAYBE_UNUSED(timebuf);

    v.output_slice(string("analytical_cube_slice") + "_" + string(timebuf) + ".dat",
                   'Y', 0.05, 'X', 0.0, 0.1, 0.005, 'Z', 0.0, 0.1, 0.005);

    cout << v.solution(point(0.01,0.01,0.01)) << endl;

    size_t maxi = 3;
    double h = 0.04, x = 0.01, y = 0.01, z = 0.01;
    double diff_a = 0.0, diff_ab = 0.0;
    for(size_t i = 0; i < maxi; i++)
    {
        y = 0.01;
        for(size_t j = 0; j < maxi; j++)
        {
            z = 0.01;
            for(size_t k = 0; k < maxi; k++)
            {
                point p(x, y, z);
                cvector3 a = v.solution(p);
                cvector3 b = func_true(p);
                cvector3 c = a - b;
                cvector3 d = a;
                diff_ab += abs(c.x * c.x + c.y * c.y + c.z * c.z);
                diff_a += abs(d.x * d.x + d.y * d.y + d.z * d.z);
                z += h;
            }
            y += h;
        }
        x += h;
    }
    v.calculate_diff();
    cout << "Diff (C3): \t" << sqrt(diff_ab / diff_a) << endl;
}

#endif // ANALYTICAL_CUBE
