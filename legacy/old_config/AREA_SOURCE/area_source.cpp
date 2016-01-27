#include "problems.h"

#if defined AREA_SOURCE
#if defined VFEM_USE_NONHOMOGENEOUS_FIRST || defined VFEM_USE_ANALYTICAL || defined VFEM_USE_PML
#error "Please, reconfigure!"
#endif

double SLAE_MAIN_EPSILON = 1e-7;

cvector3 func_rp(const point & p, const phys_area & phys)
{
    MAYBE_UNUSED(p);
    MAYBE_UNUSED(phys);

    return cvector3(0.0, 0.0, 0.0);
}

//string mesh_filename = "data/source/hg.msh";
//string phys_filename = "data/source/phys.txt";
string mesh_filename = "data/source/test.msh";
string phys_filename = "data/source/test.txt";
string tecplot_filename = "source.plt";

void postprocessing(VFEM & v, char * timebuf)
{
    MAYBE_UNUSED(v);
    MAYBE_UNUSED(timebuf);

    v.output_slice(string("source_xy") + "_" + string(timebuf) + ".dat",
                   'Z', 0.0, 'X', -100, 100, 4.0, 'Y', -100, 100, 4.0);
    v.output_slice(string("source_xz") + "_" + string(timebuf) + ".dat",
                   'Y', 0.0, 'X', -100, 100, 4.0, 'Z', -100, 100, 4.0);
}

#endif // AREA_SOURCE
