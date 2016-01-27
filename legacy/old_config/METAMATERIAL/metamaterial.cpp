#include "problems.h"

#if defined METAMATERIAL
#if !defined VFEM_USE_NONHOMOGENEOUS_FIRST || defined VFEM_USE_ANALYTICAL || defined VFEM_USE_PML
#error "Please, reconfigure!"
#endif

double SLAE_SURF_EPSILON = 1e-15;
double SLAE_MAIN_EPSILON = 1e-14;

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
    case 176:
        return cvector3(0.0, 0.0, -1.0);
    case 177:
        return cvector3(0.0, 1.0, 0.0);
    case 178:
        return cvector3(0.0, 0.0, 1.0);
    case 179:
        return cvector3(0.0, -1.0, 0.0);
    }
    cerr << "Unknown bound!" << endl;
    return cvector3(0.0, 0.0, 0.0);
}

string mesh_filename = "data/metamaterial/crystal_2dim_v2.msh";
string phys_filename = "data/metamaterial/crystal_2dim_v2.txt";
string tecplot_filename = "crystal_2dim_v2.plt";

void postprocessing(VFEM & v, char * timebuf)
{
    MAYBE_UNUSED(v);
    MAYBE_UNUSED(timebuf);
}

#endif // METAMATERIAL
