#include "triangle.h"

// ============================================================================

// Индексы для построения базисных функций на треугольниках
namespace tr_basis_indexes
{
    // Edges (Ребра) // k, l : k < l
    const size_t ind_e[3][2] =
    {
        { 0, 1 },
        { 0, 2 },
        { 1, 2 }
    };
}

// ============================================================================

matrix_t<double> triangle_base::M() const
{
    assert(false);
    return matrix_t<double>(1, 1);
}

array_t<complex<double> > triangle_base::rp(eval_func, void *)
{
    assert(false);
    return array_t<complex<double> >(1);
}

// ============================================================================

triangle_full::triangle_full(const triangle_base & other) : triangle_basic(other)
{
    basis = NULL;
}

void triangle_full::init(const basis_type * basis)
{
    init_3d(basis->tr_int);
    this->basis = basis;
}

vector3 triangle_full::w(size_t i, const point & p) const
{
    using namespace tr_basis_indexes;
    assert(i < basis->tr_bf_num);

    point p_loc = to_local(p);

    // Первый неполный
    if(i < 3)
    {
        return lambda(ind_e[i][0], p_loc) * grad_lambda(ind_e[i][1]) -
               lambda(ind_e[i][1], p_loc) * grad_lambda(ind_e[i][0]);
    }
    // Первый полный
    else if(i < 6)
    {
        size_t ii = i - 3;
        return lambda(ind_e[ii][0], p_loc) * grad_lambda(ind_e[ii][1]) +
               lambda(ind_e[ii][1], p_loc) * grad_lambda(ind_e[ii][0]);
    }
    // Второй неполный
    else if(i < 7)
    {
        return lambda(1, p_loc) * lambda(2, p_loc) * grad_lambda(0) +
               lambda(0, p_loc) * lambda(2, p_loc) * grad_lambda(1) -
               2.0 * lambda(0, p_loc) * lambda(1, p_loc) * grad_lambda(2);
    }
    else if(i < 8)
    {
        return lambda(1, p_loc) * lambda(2, p_loc) * grad_lambda(0) -
               2.0 * lambda(0, p_loc) * lambda(2, p_loc) * grad_lambda(1) +
               lambda(0, p_loc) * lambda(1, p_loc) * grad_lambda(2);
    }
    // Второй полный
    else if(i < 9)
    {
        return lambda(1, p_loc) * lambda(2, p_loc) * grad_lambda(0) +
               lambda(0, p_loc) * lambda(2, p_loc) * grad_lambda(1) +
               lambda(0, p_loc) * lambda(1, p_loc) * grad_lambda(2);
    }
    else if(i < 12)
    {
        size_t ii = i - 9;
        return lambda(ind_e[ii][1], p_loc) * (2.0 * lambda(ind_e[ii][0], p_loc) - lambda(ind_e[ii][1], p_loc)) * grad_lambda(ind_e[ii][0]) -
               lambda(ind_e[ii][0], p_loc) * (2.0 * lambda(ind_e[ii][1], p_loc) - lambda(ind_e[ii][0], p_loc)) * grad_lambda(ind_e[ii][1]);
    }

    return vector3();
}

double triangle_full::integrate_w(size_t i, size_t j) const
{
    const triangle_integration * tr_integration = &(basis->tr_int);
    double result = 0.0;
    for(size_t k = 0; k < tr_integration->get_gauss_num(); k++)
        result += tr_integration->get_gauss_weight(k) *
                  w(i, get_gauss_point(k)) *
                  w(j, get_gauss_point(k));
    return result * get_jacobian();
}

complex<double> triangle_full::integrate_fw(eval_func func, size_t i, void * data)
{
    const triangle_integration * tr_integration = &(basis->tr_int);
    complex<double> result(0.0, 0.0);
    for(size_t k = 0; k < tr_integration->get_gauss_num(); k++)
        result += tr_integration->get_gauss_weight(k) *
                  func(get_gauss_point(k), get_phys_area(), data) *
                  w(i, get_gauss_point(k));
    return result * get_jacobian();
}

matrix_t<double> triangle_full::M() const
{
    matrix_t<double> matr(basis->tr_bf_num, basis->tr_bf_num);
    for(size_t i = 0; i < basis->tr_bf_num; i++)
        for(size_t j = 0; j <= i; j++)
            matr[j][i] = matr[i][j] = integrate_w(i, j);
    return matr;
}

array_t<complex<double> >
triangle_full::rp(eval_func func, void * data)
{
    array_t<complex<double> > arr(basis->tr_bf_num);
    for(size_t i = 0; i < basis->tr_bf_num; i++)
        arr[i] = integrate_fw(func, i, data);
    return arr;
}

// ============================================================================
