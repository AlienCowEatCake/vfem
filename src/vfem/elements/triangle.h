#if !defined(TRIANGLE_H_INCLUDED)
#define TRIANGLE_H_INCLUDED

#include "../common/common.h"
#include "../common/config.h"
#include "../elements/edge.h"
#include "../vfem/phys.h"

typedef cvector3(* eval_func)(const point &, phys_area &, void *);

// Класс треугольник (обычный)
class triangle_base : public virtual triangle_basic<point, edge, face, phys_area>
{
public:
    virtual void init(const basis_type *) {}
    virtual matrix_t<double> M() const;
    virtual array_t<complex<double> > rp(eval_func, void *);
};

// Класс треугольник (полный, для работы с первыми неоднородными краевыми)
class triangle_full : public triangle_base, private triangle_basic_3d<edge, face, phys_area>
{
public:
    triangle_full(const triangle_base & other = triangle_base());

    virtual void init(const basis_type * basis);
    // Локальная матрица массы
    virtual matrix_t<double> M() const;
    // Локальная правая часть
    virtual array_t<complex<double> > rp(eval_func func, void * data);

protected:

    // Параметры базиса
    const basis_type * basis;
    // Базисные функции
    vector3 w(size_t i, const point & p) const;

    double integrate_w(size_t i, size_t j) const;
    complex<double> integrate_fw(eval_func func, size_t i, void * data);
};

typedef triangle_base triangle;

#endif // TRIANGLE_H_INCLUDED
