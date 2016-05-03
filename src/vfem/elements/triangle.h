#if !defined(TRIANGLE_H_INCLUDED)
#define TRIANGLE_H_INCLUDED

#include "../common/common.h"
#include "../common/cubatures.h"
#include "../config/config.h"
#include "../elements/edge.h"
#include "../vfem/phys.h"

typedef cvector3(* eval_func)(const point &, const phys_area &, void *);

// Индексы для построения базисных функций на треугольниках
namespace tr_basis_indexes
{
    // Edges (Ребра) // k, l : k < l
    extern const size_t ind_e[3][2];
}

// Класс треугольник (обычный)
class triangle_base
{
public:
    triangle_base();

    point * nodes[3];   // Узлы
    edge * edges[3];    // Ребра
    const point & get_node(size_t i) const;
    const edge & get_edge(size_t i) const;
    face * faces;       // Грани
    const face & get_face() const;

    phys_area * phys;   // Физическая область
    const phys_area & get_phys_area() const;

    virtual void init(const basis_type *) {}
    virtual matrix_t<double> M() const;
    virtual array_t<complex<double> > rp(eval_func, void *) const;
};

// Класс треугольник (полный, для работы с первыми неоднородными краевыми)
class triangle_full : public triangle_base
{
public:
    triangle_full(const triangle_base & other = triangle_base());

    virtual void init(const basis_type * basis);
    // Локальная матрица массы
    virtual matrix_t<double> M() const;
    // Локальная правая часть
    virtual array_t<complex<double> > rp(eval_func func, void * data) const;

protected:
    // Матрица L-координат
    matrix_t<double> L;
    // L-координаты
    double lambda(size_t i, const point & p) const;
    // Градиент L-координат в глобальных координатах
    vector3 grad_lambda(size_t i) const;

    // Параметры базиса
    const basis_type * basis;
    // Базисные функции
    vector3 w(size_t i, const point & p) const;

    // Матрица перехода между локальной и глобальной с.к.
    matrix_t<double, 3, 3> transition_matrix;
    point to_local(const point & p) const;
    point to_global(const point & p) const;
    vector3 to_global(const vector3 & v) const;

    // Точки Гаусса
    array_t<point> gauss_points;
    // Якобиан
    double jacobian;

    double integrate_w(size_t i, size_t j) const;
    complex<double> integrate_fw(eval_func func, size_t i, void * data) const;
};

typedef triangle_base triangle;

#endif // TRIANGLE_H_INCLUDED
