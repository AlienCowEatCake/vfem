#include "problems.h"

#if defined AREA_2LAYERS_LOOP_UNIVERSAL_PML
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
    if(fe->phys->gmsh_num >= 20)
        return true;
    return false;
}

class pml_config
{
public:
    double begin;
    double width;
    complex<double> chi;
    double m;
    pml_config()
    {
        begin = 600.0;
        width = 100.0;
        chi.real(4.0);
        chi.imag(1.0);
        m = 3.0;
        load();
    }
    void load()
    {
        ifstream ifs;
        string name;
        double tmp, tmp2;

        name = "pml_chi.txt";
        ifs.open(name.c_str(), ios::in);
        ifs >> tmp >> tmp2;
        if(!ifs.good())
        {
            cerr << "Error in " << __FILE__ << ":" << __LINE__
                 << " while reading file " << name << endl;
            //throw IO_FILE_ERROR;
        }
        else
        {
            chi.real(tmp);
            chi.imag(tmp2);
        }
        ifs.close();

        name = "pml_m.txt";
        ifs.open(name.c_str(), ios::in);
        ifs >> tmp;
        if(!ifs.good())
        {
            cerr << "Error in " << __FILE__ << ":" << __LINE__
                 << " while reading file " << name << endl;
            //throw IO_FILE_ERROR;
        }
        else
            m = tmp;
        ifs.close();

        name = "pml_begin.txt";
        ifs.open(name.c_str(), ios::in);
        ifs >> tmp;
        if(!ifs.good())
        {
            cerr << "Error in " << __FILE__ << ":" << __LINE__
                 << " while reading file " << name << endl;
            //throw IO_FILE_ERROR;
        }
        else
            begin = tmp;
        ifs.close();

        name = "pml_width.txt";
        ifs.open(name.c_str(), ios::in);
        ifs >> tmp;
        if(!ifs.good())
        {
            cerr << "Error in " << __FILE__ << ":" << __LINE__
                 << " while reading file " << name << endl;
            //throw IO_FILE_ERROR;
        }
        else
            width = tmp;
        ifs.close();
    }
};

// Все данные будем читать из конфига
static pml_config config;

cvector3 get_s(const point & p, const finite_element * fe, const phys_pml_area * phys_pml)
{
    if(!is_pml(p, fe))
        return cvector3(1.0, 1.0, 1.0);

    /// Нефиг пост-PML растягивать!
    if(fabs(p.x) > config.begin + config.width + 1.0 || fabs(p.y) > config.begin + config.width + 1.0 || fabs(p.z) > config.begin + config.width + 1.0)
        return cvector3(1.0, 1.0, 1.0);

    double m = config.m;
    complex<double> chi = config.chi;

    double pml_thickness_x = config.width;
    double pml_thickness_y = config.width;
    double pml_thickness_z = config.width;

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

string tecplot_filename = "area_2layers_loop_universal_pml.plt";

string phys_filename_pml = "data/area_2layers_loop_pml/std-1.txt";
string phys_filename_nonpml = "data/area_2layers_loop_pml/std-2.txt";
#if !defined SMALL_MESH
string mesh_filename = "data/area_2layers_loop_universal_pml/universal_full.msh";
#else
string mesh_filename = "data/area_2layers_loop_universal_pml/universal_small.msh";
#endif
#if defined VFEM_USE_PML
string phys_filename = phys_filename_pml;
#else
string phys_filename = phys_filename_nonpml;
#endif
string slae_dump_filename = "area_2layers_loop_universal_pml_slae.txt";

void postprocessing(VFEM & v, char * timebuf)
{
    MAYBE_UNUSED(v);
    MAYBE_UNUSED(timebuf);
#if !defined SMALL_MESH
    v.output_slice(string("area_2layers_loop_universal_pml") + "_" + string(timebuf) + ".dat",
                   'Y', 0.0, 'X', -700, 700, 20.0, 'Z', -700, 700, 20.0);

    v.output_slice(string("area_2layers_loop_universal_pml_z=10") + "_" + string(timebuf) + ".dat",
                   'Z', 10.0, 'X', -500, 500, 10.0, 'Y', -500, 500, 10.0);
    v.output_slice(string("area_2layers_loop_universal_pml_z=-10") + "_" + string(timebuf) + ".dat",
                   'Z', -10.0, 'X', -500, 500, 10.0, 'Y', -500, 500, 10.0);

    double y = 0, z0 = 10, z1 = -10;
    size_t n = 65000;
    double x0 = -500.0, x1 = 500.0;
    double hx = (x1 - x0) / (double)n;
#if defined VFEM_USE_PML
    ofstream ff((string("line_pml") + "_" + string(timebuf) + ".txt").c_str());
#else
    ofstream ff((string("line_std") + "_" + string(timebuf) + ".txt").c_str());
#endif
    for(size_t i = 0; i <= n; i++)
    {
        double x = x0 + (double)i * hx;
        cvector3 sol0 = v.solution(point(x, y, z0));
        cvector3 sol1 = v.solution(point(x, y, z1));
        ff << x << " " << sol0.x.real() << " " << sol0.x.imag()
                << " " << sol0.y.real() << " " << sol0.y.imag()
                << " " << sol0.z.real() << " " << sol0.z.imag()
                << " " << sol1.x.real() << " " << sol1.x.imag()
                << " " << sol1.y.real() << " " << sol1.y.imag()
                << " " << sol1.z.real() << " " << sol1.z.imag() << endl;
    }
    ff.close();

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
        if(fabs(v.fes[k].barycenter.x) <= 500 && fabs(v.fes[k].barycenter.x) >= 180 &&
           fabs(v.fes[k].barycenter.y) <= 500 && fabs(v.fes[k].barycenter.y) >= 180 &&
           fabs(v.fes[k].barycenter.z) <= 500 && fabs(v.fes[k].barycenter.z) >= 180)
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

#endif // AREA_2LAYERS_LOOP_UNIVERSAL_PML
