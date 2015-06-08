#include "problems.h"

#if defined FIELD_CUBE
#if defined VFEM_USE_ANALYTICAL
#error "Please, reconfigure!"
#endif

double SLAE_SURF_EPSILON = 1e-15;
double SLAE_MAIN_EPSILON = 1e-15;

cvector3 func_rp(const point & p, const phys_area & phys)
{
    MAYBE_UNUSED(p);
    MAYBE_UNUSED(phys);

    return cvector3(0.0, 0.0, 0.0);
}

cvector3 func_b1(const point & p, const triangle * tr)
{
    MAYBE_UNUSED(p);
    MAYBE_UNUSED(tr);

    size_t phys_num = tr->phys->gmsh_num;
    switch(phys_num)
    {
    case 30:
        return cvector3(0.0, 1.0, 0.0);
    case 29:
        return cvector3(0.0, -1.0, 0.0);
    case 32:
        return cvector3(-1.0, 0.0, 0.0);
    case 33:
        return cvector3(1.0, 0.0, 0.0);
    }
    cerr << "Unknown bound!" << endl;
    return cvector3(0.0, 0.0, 0.0);
}

string mesh_filename = "data/analytical_cube/cube_x4.msh";
string phys_filename = "data/analytical_cube/field_phys.txt";
string tecplot_filename = "field_cube.plt";

void postprocessing(VFEM & v, char * timebuf)
{
    MAYBE_UNUSED(v);
    MAYBE_UNUSED(timebuf);

    v.output_slice(string("field_cube_slice") + "_" + string(timebuf) + ".dat",
                   'Z', 0.5, 'X', 0.0, 1.0 + 1e-10, 0.01, 'Y', 0.0, 1.0 + 1e-10, 0.01);
}

#endif // FIELD_CUBE
