#include "problems.h"

#if defined LOOP_PML
#if defined VFEM_USE_NONHOMOGENEOUS_FIRST || defined VFEM_USE_ANALYTICAL
#error "Please, reconfigure!"
#endif

//#define SMALL_MESH
double SLAE_MAIN_EPSILON = 1e-10;

cvector3 func_rp(const point & p, const phys_area & phys)
{
    MAYBE_UNUSED(p);
    MAYBE_UNUSED(phys);
    return cvector3(0.0, 0.0, 0.0);
}

bool is_pml(const point & p, const finite_element * fe)
{
    MAYBE_UNUSED(p);
    if(fe->phys->gmsh_num != 2072 && fe->phys->gmsh_num != 2076)
        return false;
    return true;
}

cvector3 get_s(const point & p, const finite_element * fe, const phys_pml_area * phys_pml)
{
    if(!is_pml(p, fe))
        return cvector3(1.0, 1.0, 1.0);

    double m = 3;
//    complex<double> chi(3.7, 2.7);

    static complex<double> chi(-1, -1);
    if(chi.real() < 0)
    {
        ifstream chi_st;
        char chi_name[] = "chi.txt";
        chi_st.open(chi_name, ios::in);
        if(!chi_st.good())
        {
            cerr << "Error in " << __FILE__ << ":" << __LINE__
                 << " while reading file " << chi_name << endl;
            throw IO_FILE_ERROR;
        }
        double chi_re, chi_im;
        chi_st >> chi_re >> chi_im;
        chi = complex<double>(chi_re, chi_im);
        chi_st.close();
    }

    double pml_thickness_x = 100.0;
    double pml_thickness_y = 100.0;
    double pml_thickness_z = 100.0;

    double pml_distance_x = 0;
    double pml_distance_y = 0;
    double pml_distance_z = 0;
    if(p.x > phys_pml->x1) pml_distance_x = fabs(p.x - phys_pml->x1);
    if(p.x < phys_pml->x0) pml_distance_x = fabs(p.x - phys_pml->x0);
    if(p.y > phys_pml->y1) pml_distance_y = fabs(p.y - phys_pml->y1);
    if(p.y < phys_pml->y0) pml_distance_y = fabs(p.y - phys_pml->y0);
    if(p.z > phys_pml->z1) pml_distance_z = fabs(p.z - phys_pml->z1);
    if(p.z < phys_pml->z0) pml_distance_z = fabs(p.z - phys_pml->z0);

    // Коэффициенты для проекции
    double dist = sqrt(pml_distance_x * pml_distance_x + pml_distance_y * pml_distance_y + pml_distance_z * pml_distance_z);
    double cx = pml_distance_x / dist;
    double cy = pml_distance_y / dist;
    double cz = pml_distance_z / dist;
    if(is_fpu_error(cx)) cx = 0;
    if(is_fpu_error(cy)) cy = 0;
    if(is_fpu_error(cz)) cz = 0;

    // Уравнение прямой
    // (x - x1) / (x2 - x1) = (y - y1) / (y2 - y1) = (z - z1) / (z2 - z1)
    // Ищем по какой координате больше прошло, та и будет браться как максимальная
    double tx = 0, ty = 0, tz = 0;
    if(pml_distance_x > pml_distance_y && pml_distance_x > pml_distance_z) // x
    {
        tx = pml_thickness_x;
        ty = tx / pml_distance_x * pml_distance_y;
        tz = tx / pml_distance_x * pml_distance_z;
    }
    else if(fabs(pml_distance_y) > fabs(pml_distance_x) && fabs(pml_distance_y) > fabs(pml_distance_z)) // y
    {
        ty = pml_thickness_y;
        tx = ty / pml_distance_y * pml_distance_x;
        tz = ty / pml_distance_y * pml_distance_z;
    }
    else    // z
    {
        tz = pml_thickness_z;
        tx = tz / pml_distance_z * pml_distance_x;
        ty = tz / pml_distance_z * pml_distance_y;
    }
    if(is_fpu_error(tx)) tx = 0;
    if(is_fpu_error(ty)) ty = 0;
    if(is_fpu_error(tz)) tz = 0;
    double tick = sqrt(tx * tx + ty * ty + tz * tz);

    double power = pow(dist / tick, m);
    return cvector3(1.0 + chi * power * cx, 1.0 + chi * power * cy, 1.0 + chi * power * cz);
}

//string phys_filename = "data/loop_pml/1.txt";
//string mesh_filename = "data/loop_pml/1.msh";
string tecplot_filename = "loop_pml.plt";

string phys_filename_pml = "data/loop_pml/2-1.txt";
string phys_filename_nonpml = "data/loop_pml/2-2.txt";
#if !defined SMALL_MESH
string mesh_filename = "data/loop_pml/3.msh";
#else
string mesh_filename = "data/loop_pml/3-small.msh";
#endif
#if defined VFEM_USE_PML
string phys_filename = phys_filename_pml;
#else
string phys_filename = phys_filename_nonpml;
#endif
string slae_dump_filename = "loop_pml_slae.txt";

void postprocessing(VFEM & v, char * timebuf)
{
    MAYBE_UNUSED(v);
    MAYBE_UNUSED(timebuf);
#if !defined SMALL_MESH
    v.output_slice(string("loop_pml_slice") + "_" + string(timebuf) + ".dat",
                   'Z', 0.0, 'X', -700, 700, 70, 'Y', -700, 700, 70);
#if !defined VFEM_USE_PML
    v.slae.dump_x(slae_dump_filename);
#else
    complex<double> * anal = new complex<double> [v.slae.n];
    ifstream a;
    a.open(slae_dump_filename.c_str(), ios::in);
    for(size_t i = 0; i < v.slae.n; i++)
        a >> anal[i];
    a.close();

    double diff = 0.0, norm = 0.0;
    for(size_t k = 0; k < v.fes.size(); k++)
    {
        if(fabs(v.fes[k].barycenter.x) <= 600 && fabs(v.fes[k].barycenter.x) >= 180 &&
           fabs(v.fes[k].barycenter.y) <= 600 && fabs(v.fes[k].barycenter.y) >= 180 &&
           fabs(v.fes[k].barycenter.z) <= 600 && fabs(v.fes[k].barycenter.z) >= 180)
        {
            array_t<complex<double>, basis::tet_bf_num> q_loc, q_loc_true;
            for(size_t i = 0; i < basis::tet_bf_num; i++)
            {
                size_t dof = v.fes[k].dof[i];
                q_loc[i] = v.slae.x[dof];
                q_loc_true[i] = anal[dof];
            }
            diff += v.fes[k].diff_normL2(q_loc, q_loc_true);
            norm += v.fes[k].normL2(q_loc_true);
        }
    }
    cout << "Diff (L2): \t" << sqrt(diff / norm) << endl;
    delete [] anal;
#endif
#endif
}

#endif // LOOP_PML
