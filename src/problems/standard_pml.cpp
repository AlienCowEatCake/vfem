#include "problems.h"

#if defined VFEM_USE_PML

// Проверка, PML или нет
bool is_pml(const point & p, const finite_element * fe, const phys_pml_area * phys_pml)
{
    MAYBE_UNUSED(p);
    if(phys_pml->params.find(fe->get_phys_area().gmsh_num) != phys_pml->params.end())
        return true;
    return false;
}

// Коэффициент S для PML
cvector3 get_s(const point & p, const finite_element * fe, const phys_pml_area * phys_pml)
{
    if(!is_pml(p, fe, phys_pml))
        return cvector3(1.0, 1.0, 1.0);

    const pml_config_parameter * conf = & (phys_pml->params.find(fe->get_phys_area().gmsh_num)->second);
    double m = conf->m, width = conf->width;
    complex<double> chi = conf->chi;

    /// Нефиг пост-PML растягивать!
    if(p.x > phys_pml->x1 + width + 1.0 || p.x < phys_pml->x0 - width - 1.0 || p.y > phys_pml->y1 + width + 1.0 ||
       p.y < phys_pml->y0 - width - 1.0 || p.z > phys_pml->z1 + width + 1.0 || p.z < phys_pml->z0 - width - 1.0)
        return cvector3(1.0, 1.0, 1.0);

    double pml_thickness_x = width;
    double pml_thickness_y = width;
    double pml_thickness_z = width;

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

#endif
