#include "problems.h"

#if defined AREA_POINT_SOURCE
#if defined VFEM_USE_NONHOMOGENEOUS_FIRST || defined VFEM_USE_ANALYTICAL || defined VFEM_USE_PML
#error "Please, reconfigure!"
#endif

double SLAE_MAIN_EPSILON = 1e-5;

cvector3 func_rp(const point & p, const phys_area & phys)
{
    MAYBE_UNUSED(p);
    MAYBE_UNUSED(phys);

    return cvector3(0.0, 0.0, 0.0);
}

string mesh_filename = "data/area_point_source/area_symm.msh";
string phys_filename = "data/area_point_source/phys.txt";
string tecplot_filename = "area_point_source.plt";

void postprocessing(VFEM & v, char * timebuf)
{
    MAYBE_UNUSED(v);
    MAYBE_UNUSED(timebuf);

    v.output_slice(string("area_point_source_slice") + "_" + string(timebuf) + ".dat",
                   'Y', 5.0, 'X', -200.0, 200.0, 8.0, 'Z', -350.0, 0.0, 7.0);
}

#endif // AREA_POINT_SOURCE
