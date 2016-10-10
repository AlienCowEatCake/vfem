#include "phys.h"

namespace consts
{
    const double c = 299792458.0;                // Скорость света
    const double mu0 = 4.0e-7 * M_PI;            // Магн. пр. вакуума
    const double epsilon0 = 1.0 / (mu0 * c * c); // Диэл. пр. вакуума
}

phys_area::phys_area()
    : omega(0.0)
    , mu(0.0)
    , epsilon(0.0)
    , gmsh_num(0)
    , type_of_elem(0)
    , type_of_bounds(0)
    , J0(0.0)
    , E0(0.0)
{
    for(size_t i = 0; i < sigma.size(); i++)
        sigma[i].parse("0.0");
}

phys_id::phys_id(size_t type_of_element, size_t gmsh_num)
    : type_of_element(type_of_element)
    , gmsh_num(gmsh_num)
{}

bool phys_id::operator < (const phys_id & other) const
{
    if(type_of_element < other.type_of_element) return true;
    if(type_of_element > other.type_of_element) return false;
    return gmsh_num < other.gmsh_num;
}

pml_config_parameter::pml_config_parameter(complex<double> n_chi, double n_width, double n_m)
{
    chi = n_chi;
    width = n_width;
    m = n_m;
}
