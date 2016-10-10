#if !defined(PHYS_H_INCLUDED)
#define PHYS_H_INCLUDED

#include "../common/common.h"

namespace consts
{
    extern const double c;          // Скорость света
    extern const double mu0;        // Магн. пр. вакуума
    extern const double epsilon0;   // Диэл. пр. вакуума
}

// Структура физическая область
struct phys_area
{
    double omega;                // Циклическая частота
    double mu;                   // Магнитная проницаемость (относительная)
    array_t<evaluator_xyz<double>, 3> sigma;    // Электрическая проводимость
    double epsilon;              // Диэлектрическая проницаемость (относительная)
    size_t gmsh_num;             // Номер области в Gmsh
    size_t type_of_elem;         // Тип элементов
    size_t type_of_bounds;       // Тип краевого условия
    double J0;                   // Мощность источника
    double E0;                   // Электрическое поле от электрода
    phys_area();                 // Конструктор по умолчанию
};

// Структура для индексации физических областей
struct phys_id
{
    size_t type_of_element;
    size_t gmsh_num;
    phys_id(size_t type_of_element = 0, size_t gmsh_num = 0);
    bool operator < (const phys_id & other) const;
};

// Структура для хранения параметров PML для одной физической области
struct pml_config_parameter
{
    complex<double> chi;
    double width;
    double m;
    pml_config_parameter(complex<double> n_chi = 1.0, double n_width = 100.0, double n_m = 3.0);
};

// Структура для хранения параметров PML
struct phys_pml_area
{
    double x0, x1;
    double y0, y1;
    double z0, z1;
    map<size_t, pml_config_parameter> params;
};

#endif // PHYS_H_INCLUDED
