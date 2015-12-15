#include "problems.h"

#if defined VFEM_USE_PML

class pml_config_parameter
{
public:
    complex<double> chi;
    double width;
    double m;
    pml_config_parameter(complex<double> n_chi = 1.0, double n_width = 100.0, double n_m = 3.0);
};

class pml_config
{
public:
    map<size_t, pml_config_parameter> params;
    pml_config();
};


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

// ============================================================================

string trim(const string & str);
string to_lowercase(const string & str);

pml_config_parameter::pml_config_parameter(complex<double> n_chi, double n_width, double n_m)
{
    chi = n_chi;
    width = n_width;
    m = n_m;
}

pml_config::pml_config()
{
    string filename = "config_pml.ini";
    ifstream ifs(filename.c_str());
    if(!ifs.good())
    {
        cerr << "Error in " << __FILE__ << ":" << __LINE__
             << " while reading file " << filename << endl;
        exit(1);
    }
    ifs.close();
}

#endif
