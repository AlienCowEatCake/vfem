#include "point.h"

// Проверка, лежит ли точка внутри параллелепипеда (для дерева)
template<>
bool point_t<double>::inside(double x0, double x1, double y0, double y1, double z0, double z1) const
{
    if(x >= x0 && x <= x1 && y >= y0 && y <= y1 && z >= z0 && z <= z1)
        return true;
    return false;
}
