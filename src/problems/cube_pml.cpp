#include "problems.h"

#if defined CUBE_PML
#if defined VFEM_USE_NONHOMOGENEOUS_FIRST || (!defined VFEM_USE_PML && defined VFEM_USE_ANALYTICAL)
#error "Please, reconfigure!"
#endif

double SLAE_MAIN_EPSILON = 1e-7;

cvector3 func_rp(const point & p, const phys_area & phys)
{
    MAYBE_UNUSED(p);
    MAYBE_UNUSED(phys);
    return cvector3(0.0, 0.0, 0.0);
}

bool is_pml(const point & p, const finite_element * fe)
{
    //return false;
    MAYBE_UNUSED(p);
    if(fe->phys->gmsh_num != 106 && fe->phys->gmsh_num != 107)
        return false;
    return true;
}

cvector3 get_s(const point & p, const finite_element * fe, const phys_pml_area * phys_pml)
{
    if(!is_pml(p, fe))
        return cvector3(1.0, 1.0, 1.0);

    double m = 3;
    complex<double> chi(5, 5);
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

string phys_filename_pml = "data/cube_pml/nonpml-phys-first.txt";
string phys_filename_nonpml = "data/cube_pml/nonpml-phys-second.txt";
string slae_dump_filename = "cube_nonpml_slae.txt";

//string mesh_filename = "data/cube_pml/pml-test.msh";
string mesh_filename = "data/cube_pml/nonpml-test.msh";
#if defined VFEM_USE_PML
string phys_filename = phys_filename_pml;
#else
string phys_filename = phys_filename_nonpml;
#endif
string tecplot_filename = "cube_pml.plt";

#if defined VFEM_USE_ANALYTICAL
cvector3 func_true(const point & p)
{
    static VFEM vfem_anal;
    if(vfem_anal.fes_num == 0)
    {
        vfem_anal.input_phys(phys_filename_nonpml);
        vfem_anal.input_mesh(mesh_filename);
        vfem_anal.slae.restore(slae_dump_filename);
    }
    return vfem_anal.solution(p);
}
#endif

void postprocessing(VFEM & v, char * timebuf)
{
    MAYBE_UNUSED(v);
    MAYBE_UNUSED(timebuf);
    v.output_slice(string("cube_pml_slice") + "_" + string(timebuf) + ".dat",
                   'Z', 0.0, 'X', -700, 700, 20.0, 'Y', -700, 700, 20.0);
#if !defined VFEM_USE_PML
    v.slae.dump(slae_dump_filename);
#elif defined VFEM_USE_ANALYTICAL
    double diff = 0.0, norm = 0.0;
    for(size_t k = 0; k < v.fes_num; k++)
    {
        if(fabs(v.fes[k].barycenter.x) <= 600 && fabs(v.fes[k].barycenter.x) >= 10 &&
           fabs(v.fes[k].barycenter.y) <= 600 && fabs(v.fes[k].barycenter.y) >= 10 &&
           fabs(v.fes[k].barycenter.z) <= 600 && fabs(v.fes[k].barycenter.z) >= 10)
        {
            array_t<complex<double>, basis::tet_bf_num> q_loc;
            for(size_t i = 0; i < basis::tet_bf_num; i++)
            {
                size_t dof = v.fes[k].dof[i];
                q_loc[i] = v.slae.x[dof];
            }
            diff += v.fes[k].diff_normL2(q_loc, func_true);
            norm += v.fes[k].normL2(func_true);
        }
    }
    cout << "Diff (L2): \t" << sqrt(diff / norm) << endl;
#endif
}

#endif // CUBE_PML
