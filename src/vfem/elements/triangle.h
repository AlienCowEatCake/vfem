#if !defined(TRIANGLE_H_INCLUDED)
#define TRIANGLE_H_INCLUDED

#include "../common/common.h"
#include "../common/config.h"
#include "../elements/edge.h"
#include "../vfem/phys.h"

typedef cvector3(* eval_func)(const point &, phys_area &, void *);

// Класс треугольник (обычный)
struct triangle_base : public triangle_basic<point, edge, face, phys_area>
{
    virtual void init(const basis_type *) {}
    virtual matrix_t<double> M() const;
    virtual array_t<complex<double> > rp(eval_func, void *);
};

// Класс треугольник (полный, для работы с первыми неоднородными краевыми)
struct triangle_full : public triangle_base
{
public:
    triangle_full(const triangle_base & other = triangle_base());

    virtual void init(const basis_type * basis);
    // Локальная матрица массы
    virtual matrix_t<double> M() const;
    // Локальная правая часть
    virtual array_t<complex<double> > rp(eval_func func, void * data);

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
    complex<double> integrate_fw(eval_func func, size_t i, void * data);
};

typedef triangle_base triangle;

#endif // TRIANGLE_H_INCLUDED
