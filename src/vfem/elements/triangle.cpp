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

triangle_full::triangle_full(const triangle_base & other) : triangle_base(other)
{
    basis = NULL;
    jacobian = 0.0;
}

void triangle_full::init(const basis_type * basis)
{
    const triangle_integration * tr_integration = &(basis->tr_int);

    // Построение локальной системы координат
    vector3 g1(get_node(0), get_node(1));
    vector3 g2(get_node(0), get_node(2));
    vector3 e1 = g1 / g1.norm();
    vector3 e2 = g2 - ((g1 * g2) / (g1 * g1)) * g1;
    e2 = e2 / e2.norm();
    vector3 e3 = g1.cross(e2);
    e3 = e3 / e3.norm();

    for(size_t i = 0; i < 3; i++)
    {
        transition_matrix[0][i] = e1[i];
        transition_matrix[1][i] = e2[i];
        transition_matrix[2][i] = e3[i];
    }

    point local_nodes[3];
    local_nodes[0] = point(0.0, 0.0, 0.0);
    local_nodes[1] = point(g1.norm(), 0.0, 0.0);
    local_nodes[2] = (transition_matrix * g2).to_point();

    matrix_t<double, 3, 3> D;
    // Формирование L-координат
    for(size_t i = 0; i < 2; i++)
        for(size_t j = 0; j < 3; j++)
            D[i][j] = local_nodes[j][i];

    for(size_t i = 0; i < 3; i++)
        D[2][i] = 1.0;

    double D_det;
    L = inverse(D, D_det);

    jacobian = fabs(D_det);

    // Точки Гаусса в локальной системе координат
    array_t<point> gauss_points_local(tr_integration->get_gauss_num());
    for(size_t j = 0 ; j < tr_integration->get_gauss_num(); j++)
    {
        for(size_t i = 0; i < 2; i++)
        {
            gauss_points_local[j][i] = 0.0;
            for(size_t k = 0; k < 3; k++)
                gauss_points_local[j][i] += D[i][k] * tr_integration->get_gauss_point_master(j, k);
        }
        gauss_points_local[j][2] = 0.0;
    }

    // Точки Гаусса в глобальной системе координат
    gauss_points.resize(tr_integration->get_gauss_num());
    for(size_t i = 0; i < tr_integration->get_gauss_num(); i++)
        gauss_points[i] = to_global(gauss_points_local[i]);

    // Инициализируем параметры базиса
    this->basis = basis;
}

point triangle_full::to_local(const point & p) const
{
    point shift = p;
    // Cдвиг
    for(size_t i = 0; i < 3; i++)
        shift[i] -= get_node(0)[i];
    point turn;
    // Поворот
    for(size_t i = 0; i < 3; i++)
    {
        turn[i] = 0.0;
        for(size_t j = 0; j < 3; j++)
            turn[i] += transition_matrix[i][j] * shift[j];
    }
    return turn;
}

point triangle_full::to_global(const point & p) const
{
    point turn;
    // Поворот
    for(size_t i = 0; i < 3; i++)
    {
        turn[i] = 0.0;
        for(size_t j = 0; j < 3; j++)
            turn[i] += transition_matrix[j][i] * p[j];
    }
    point shift;
    // Сдвиг
    for(size_t i = 0; i < 3; i++)
        shift[i] = turn[i] + get_node(0)[i];
    return shift;
}

vector3 triangle_full::to_global(const vector3 & v) const
{
    vector3 result;
    // Поворот
    for(size_t i = 0; i < 3; i++)
    {
        result[i] = 0;
        for(size_t j = 0; j < 3; j++)
            result[i] += transition_matrix[j][i] * v[j];
    }
    return result;
}

double triangle_full::lambda(size_t i, const point & p) const
{
    double result = 0.0;
    for(size_t j = 0; j < 2; j++)
        result += L[i][j] * p[j];
    result += L[i][2];
    return result;
}

vector3 triangle_full::grad_lambda(size_t i) const
{
    vector3 grad(L[i][0], L[i][1], 0.0);
    return to_global(grad);
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
                  w(i, gauss_points[k]) *
                  w(j, gauss_points[k]);
    return result * jacobian;
}

complex<double> triangle_full::integrate_fw(eval_func func, size_t i, void * data)
{
    const triangle_integration * tr_integration = &(basis->tr_int);
    complex<double> result(0.0, 0.0);
    for(size_t k = 0; k < tr_integration->get_gauss_num(); k++)
        result += tr_integration->get_gauss_weight(k) *
                  func(gauss_points[k], get_phys_area(), data) *
                  w(i, gauss_points[k]);
    return result * jacobian;
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
