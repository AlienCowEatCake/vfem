#include "problems.h"

#if defined ANALYTICAL_CUBE
#if !defined VFEM_USE_ANALYTICAL
#error "Please, reconfigure!"
#endif

double SLAE_SURF_EPSILON = 1e-15;
double SLAE_MAIN_EPSILON = 1e-15;

cvector3 func_true(const point & p)
{
    MAYBE_UNUSED(p);

    //return cvector3(p.y + p.z, p.x + p.z, p.x + p.y);

    complex<double> x(exp(- (0.5 - p.y) * (0.5 - p.y) - (0.5 - p.z) * (0.5 - p.z)));
    complex<double> y(exp(- (0.5 - p.x) * (0.5 - p.x) - (0.5 - p.z) * (0.5 - p.z)));
    complex<double> z(exp(- (0.5 - p.x) * (0.5 - p.x) - (0.5 - p.y) * (0.5 - p.y)));
    return cvector3(x, y, z);
}

cvector3 func_rp(const point & p, const phys_area & phys)
{
    MAYBE_UNUSED(p);
    MAYBE_UNUSED(phys);

    complex<double> k2(- phys.omega * phys.omega * phys.epsilon, phys.omega * phys.sigma);
    //return k2 * cvector3(p.y + p.z, p.x + p.z, p.x + p.y);

    // http://www.wolframalpha.com/input/?i=curl%28curl%28%7Bexp%28-%280.5-y%29%5E2-%280.5-z%29%5E2%29%2Cexp%28-%280.5-x%29%5E2-%280.5-z%29%5E2%29%2Cexp%28-%280.5-x%29%5E2+-%280.5-y%29%5E2%29%7D%29+%29
    return cvector3(
                exp(-0.5 + p.y - p.y * p.y + p.z - p.z * p.z) * (2.0 + 4.0 * p.y - 4.0 * p.y * p.y + 4.0 * p.z - 4.0 * p.z * p.z) / phys.mu + k2 * exp(- (0.5 - p.y) * (0.5 - p.y) - (0.5 - p.z) * (0.5 - p.z)),
                exp(-0.5 + p.x - p.x * p.x + p.z - p.z * p.z) * (2.0 + 4.0 * p.x - 4.0 * p.x * p.x + 4.0 * p.z - 4.0 * p.z * p.z) / phys.mu + k2 * exp(- (0.5 - p.x) * (0.5 - p.x) - (0.5 - p.z) * (0.5 - p.z)),
                exp(-0.5 + p.x - p.x * p.x + p.y - p.y * p.y) * (2.0 + 4.0 * p.x - 4.0 * p.x * p.x + 4.0 * p.y - 4.0 * p.y * p.y) / phys.mu + k2 * exp(- (0.5 - p.x) * (0.5 - p.x) - (0.5 - p.y) * (0.5 - p.y))
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

string mesh_filename = "data/analytical_cube/cube_x4.msh";
string phys_filename = "data/analytical_cube/phys.txt";
string tecplot_filename = "analytical_cube.plt";

void postprocessing(VFEM & v, char * timebuf)
{
    MAYBE_UNUSED(v);
    MAYBE_UNUSED(timebuf);

    v.output_slice(string("analytical_cube_slice") + "_" + string(timebuf) + ".dat",
                   'Y', 0.5, 'X', 0.0, 1.0 + 1e-10, 0.01, 'Z', 0.0, 1.0 + 1e-10, 0.01);
    v.calculate_diff();
}

#endif // ANALYTICAL_CUBE
