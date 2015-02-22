#if !defined TRIANGLE_H_INCLUDED
#define TRIANGLE_H_INCLUDED

#include "../common/common.h"
#include "../common/matrix.h"
#include "../geometry/point.h"
#include "../geometry/vector3.h"
#include "../elements/edge.h"
#include "../vfem/phys.h"

// Класс треугольник (обычный)
class triangle_base
{
public:
    triangle_base();

    point * nodes[3];   // Узлы
    edge * edges[3];    // Ребра
    const point & get_node(size_t i) const;
    const edge & get_edge(size_t i) const;

    phys_area * phys;   // Физическая область
    const phys_area & get_phys_area() const;
};

// Класс треугольник (полный, для работы с первыми неоднородными краевыми)
class triangle_full : public triangle_base
{
public:
    void init();
    matrix6 M() const;  // Локальная матрица массы
    carray6 rp(cvector3(*func)(const point &, const triangle_full *)) const; // Локальная правая часть

    size_t edges_surf[3];   // Номера ребер в массиве по границе

protected:
    matrix3 L;  // Матрица L-координат
    double lambda(size_t i, const point & p) const; // L-координаты

    vector3 w(size_t i, const point & p) const;     // Базисные функции

    matrix3 transition_matrix;  // Матрица перехода между локальной и глобальной с.к.
    point to_local(const point & p) const;
    point to_global(const point & p) const;
    vector3 to_global(const vector3 & v) const;

    point gauss_points[3];      // Точки Гаусса
    double gauss_weights[3];    // Веса
    double jacobian;            // Якобиан

    double integrate_w(size_t i, size_t j) const;
    complex<double> integrate_fw(cvector3(*func)(const point &, const triangle_full *), size_t i) const;
};

#endif // TRIANGLE_H_INCLUDED
