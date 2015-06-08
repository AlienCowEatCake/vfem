#if !defined PHYS_H_INCLUDED
#define PHYS_H_INCLUDED

#include "../common/common.h"

namespace consts
{
    static const double c = 299792458.0;                // Скорость света
    static const double mu0 = 4.0e-7 * M_PI;            // Магн. пр. вакуума
    static const double epsilon0 = 1.0 / (mu0 * c * c); // Диэл. пр. вакуума
}

// Класс физическая область
class phys_area
{
public:
    double omega;           // Циклическая частота
    double mu;              // Магнитная проницаемость (относительная)
    double sigma;           // Электрическая проводимость
    double epsilon;         // Диэлектрическая проницаемость (относительная)
    size_t gmsh_num;        // Номер области в Gmsh
    size_t type_of_elem;    // Тип элементов
    size_t type_of_bounds;  // Тип краевого условия
    phys_area()             // Конструктор по умолчанию
    {
        omega = mu = sigma = epsilon = 0.0;
        gmsh_num = type_of_elem = type_of_bounds = 0;
    }
};

// Класс для индексации физических областей
class phys_id
{
public:
    size_t type_of_element;
    size_t gmsh_num;
    phys_id()
    {
        type_of_element = 0;
        gmsh_num = 0;
    }
    phys_id(size_t type_of_element, size_t gmsh_num)
    {
        this->type_of_element = type_of_element;
        this->gmsh_num = gmsh_num;
    }
    bool operator < (const phys_id & other) const
    {
        if(type_of_element < other.type_of_element) return true;
        if(type_of_element > other.type_of_element) return false;
        return gmsh_num < other.gmsh_num;
    }
};

#endif // PHYS_H_INCLUDED
