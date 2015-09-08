#include "problems.h"

#if defined VFEM_USE_PML

bool is_pml(const point & p, const finite_element * fe)
{
    // TODO
    MAYBE_UNUSED(p);
    if(fe->phys->gmsh_num == 21 || fe->phys->gmsh_num == 22 || fe->phys->gmsh_num == 23 ||
       fe->phys->gmsh_num == 31 || fe->phys->gmsh_num == 32 || fe->phys->gmsh_num == 33)
        return true;
    return false;
}

cvector3 get_s(const point & p, const finite_element * fe, const phys_pml_area * phys_pml)
{
    // TODO
    return cvector3(1.0, 1.0, 1.0);
}

#endif
