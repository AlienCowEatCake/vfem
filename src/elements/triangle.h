#if !defined TRIANGLE_H_INCLUDED
#define TRIANGLE_H_INCLUDED

#include "../common/common.h"
#include "../common/matrix.h"
#include "../geometry/point.h"
#include "../geometry/vector3.h"
#include "../elements/edge.h"
#include "../elements/face.h"
#include "../vfem/phys.h"
#include "../elements/basis_config.h"

// Класс треугольник (обычный)
class triangle_base
{
public:
    triangle_base();

    point * nodes[3];   // Узлы
    edge * edges[3];    // Ребра
    const point & get_node(size_t i) const;
    const edge & get_edge(size_t i) const;
#if BASIS_ORDER >= 2
    face * faces;       // Грани
    const face & get_face() const;
#endif

    phys_area * phys;   // Физическая область
    const phys_area & get_phys_area() const;

    size_t dof[basis::tr_bf_num];
};

// Класс треугольник (полный, для работы с первыми неоднородными краевыми)
class triangle_full : public triangle_base
{
public:
    void init();
    // Локальная матрица массы
    matrix_t<double, basis::tr_bf_num, basis::tr_bf_num> M() const;
    // Локальная правая часть
    array_t<complex<double>, basis::tr_bf_num> rp(cvector3(*func)(const point &, const triangle_full *)) const;

    // Номера степеней свободы по границе
    size_t dof_surf[basis::tr_bf_num];

protected:
    // Матрица L-координат
    matrix_t<double, 3, 3> L;
    // L-координаты
    double lambda(size_t i, const point & p) const;
    // Градиент L-координат в глобальных координатах
    vector3 grad_lambda(size_t i) const;

    // Базисные функции
    vector3 w(size_t i, const point & p) const;

    // Матрица перехода между локальной и глобальной с.к.
    matrix_t<double, 3, 3> transition_matrix;
    point to_local(const point & p) const;
    point to_global(const point & p) const;
    vector3 to_global(const vector3 & v) const;

    // Точки Гаусса
    point gauss_points[tr_integration::gauss_num];
    // Якобиан
    double jacobian;

    double integrate_w(size_t i, size_t j) const;
    complex<double> integrate_fw(cvector3(*func)(const point &, const triangle_full *), size_t i) const;
};

#endif // TRIANGLE_H_INCLUDED
