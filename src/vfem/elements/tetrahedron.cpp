#include "tetrahedron.h"

// ============================================================================

// Индексы для построения базисных функций на тетраэдрах
namespace tet_basis_indexes
{
    // Edges (Ребра) // k, l : k < l
    const size_t ind_e[6][2] =
    {
        { 0, 1 },
        { 0, 2 },
        { 0, 3 },
        { 1, 2 },
        { 1, 3 },
        { 2, 3 }
    };
    // Faces (Грани) // j, k, l : j < k < l
    const size_t ind_f[4][3] =
    {
        { 0, 1, 2 },
        { 0, 1, 3 },
        { 0, 2, 3 },
        { 1, 2, 3 }
    };
}

// ============================================================================

vector3 tetrahedron_base::w(size_t i, const point & p) const
{
    using namespace tet_basis_indexes;
    assert(i < basis->tet_bf_num);

    // Первый неполный
    if(i < 6)
    {
        return lambda(ind_e[i][0], p) * grad_lambda(ind_e[i][1]) -
               lambda(ind_e[i][1], p) * grad_lambda(ind_e[i][0]);
    }
    // Первый полный
    else if(i < 12)
    {
        size_t ii = i - 6;
        return lambda(ind_e[ii][0], p) * grad_lambda(ind_e[ii][1]) +
               lambda(ind_e[ii][1], p) * grad_lambda(ind_e[ii][0]);
    }
    // Второй неполный
    else if(i < 16)
    {
        size_t ii = i - 12;
        return lambda(ind_f[ii][1], p) * lambda(ind_f[ii][2], p) * grad_lambda(ind_f[ii][0]) +
               lambda(ind_f[ii][0], p) * lambda(ind_f[ii][2], p) * grad_lambda(ind_f[ii][1]) -
               2.0 * lambda(ind_f[ii][0], p) * lambda(ind_f[ii][1], p) * grad_lambda(ind_f[ii][2]);
    }
    else if(i < 20)
    {
        size_t ii = i - 16;
        return lambda(ind_f[ii][1], p) * lambda(ind_f[ii][2], p) * grad_lambda(ind_f[ii][0]) -
               2.0 * lambda(ind_f[ii][0], p) * lambda(ind_f[ii][2], p) * grad_lambda(ind_f[ii][1]) +
               lambda(ind_f[ii][0], p) * lambda(ind_f[ii][1], p) * grad_lambda(ind_f[ii][2]);
    }
    // Второй полный
    else if(i < 24)
    {
        size_t ii = i - 20;
        return lambda(ind_f[ii][1], p) * lambda(ind_f[ii][2], p) * grad_lambda(ind_f[ii][0]) +
               lambda(ind_f[ii][0], p) * lambda(ind_f[ii][2], p) * grad_lambda(ind_f[ii][1]) +
               lambda(ind_f[ii][0], p) * lambda(ind_f[ii][1], p) * grad_lambda(ind_f[ii][2]);
    }
    else if(i < 30)
    {
        size_t ii = i - 24;
        return lambda(ind_e[ii][1], p) * (2.0 * lambda(ind_e[ii][0], p) - lambda(ind_e[ii][1], p)) * grad_lambda(ind_e[ii][0]) -
               lambda(ind_e[ii][0], p) * (2.0 * lambda(ind_e[ii][1], p) - lambda(ind_e[ii][0], p)) * grad_lambda(ind_e[ii][1]);
    }

    return vector3();
}

double tetrahedron_base::divw(size_t i, const point & p) const
{
    using namespace tet_basis_indexes;
    assert(i < basis->tet_bf_num);

    // Первый неполный
    if(i < 6)
    {
        return 0.0;
    }
    // Первый полный
    else if(i < 12)
    {
        size_t ii = i - 6;
        return 2.0 * grad_lambda(ind_e[ii][0]) * grad_lambda(ind_e[ii][1]);
    }
    // Второй неполный
    else if(i < 20)
    {
        return 0.0;
    }
    // Второй полный
    else if(i < 24)
    {
        size_t ii = i - 20;
        return lambda(ind_f[ii][2], p) * grad_lambda(ind_f[ii][1]) * grad_lambda(ind_f[ii][0]) +
               lambda(ind_f[ii][1], p) * grad_lambda(ind_f[ii][2]) * grad_lambda(ind_f[ii][0]) +
               lambda(ind_f[ii][2], p) * grad_lambda(ind_f[ii][0]) * grad_lambda(ind_f[ii][1]) +
               lambda(ind_f[ii][0], p) * grad_lambda(ind_f[ii][2]) * grad_lambda(ind_f[ii][1]) +
               lambda(ind_f[ii][1], p) * grad_lambda(ind_f[ii][0]) * grad_lambda(ind_f[ii][2]) +
               lambda(ind_f[ii][0], p) * grad_lambda(ind_f[ii][1]) * grad_lambda(ind_f[ii][2]);
    }
    else if(i < 30)
    {
        size_t ii = i - 24;
        return (2.0 * lambda(ind_e[ii][0], p) - lambda(ind_e[ii][1], p)) * grad_lambda(ind_e[ii][1]) * grad_lambda(ind_e[ii][0]) +
               lambda(ind_e[ii][1], p) * (2.0 * grad_lambda(ind_e[ii][0]) - grad_lambda(ind_e[ii][1])) * grad_lambda(ind_e[ii][0]) -
               (2.0 * lambda(ind_e[ii][1], p) - lambda(ind_e[ii][0], p)) * grad_lambda(ind_e[ii][0]) * grad_lambda(ind_e[ii][1]) -
               lambda(ind_e[ii][0], p) * (2.0 * grad_lambda(ind_e[ii][1]) - grad_lambda(ind_e[ii][0])) * grad_lambda(ind_e[ii][1]);
    }

    return 0.0;
}

vector3 tetrahedron_base::rotw(size_t i, const point & p) const
{
    using namespace tet_basis_indexes;
    assert(i < basis->tet_bf_num);

    // Первый неполный
    if(i < 6)
    {
        return 2.0 * grad_lambda(ind_e[i][0]).cross(grad_lambda(ind_e[i][1]));
    }
    // Первый полный
    else if(i < 12)
    {
        return vector3(0.0, 0.0, 0.0);
    }
    // Второй неполный
    else if(i < 16)
    {
        size_t ii = i - 12;
        double lambda_j = lambda(ind_f[ii][0], p);
        double lambda_k = lambda(ind_f[ii][1], p);
        double lambda_l = lambda(ind_f[ii][2], p);
        vector3 grad_lambda_j = grad_lambda(ind_f[ii][0]);
        vector3 grad_lambda_k = grad_lambda(ind_f[ii][1]);
        vector3 grad_lambda_l = grad_lambda(ind_f[ii][2]);
        return (lambda_l * grad_lambda_k + lambda_k * grad_lambda_l).cross(grad_lambda_j) +
               (lambda_l * grad_lambda_j + lambda_j * grad_lambda_l).cross(grad_lambda_k) -
               2.0 * (lambda_k * grad_lambda_j + lambda_j * grad_lambda_k).cross(grad_lambda_l);
    }
    else if(i < 20)
    {
        size_t ii = i - 16;
        double lambda_j = lambda(ind_f[ii][0], p);
        double lambda_k = lambda(ind_f[ii][1], p);
        double lambda_l = lambda(ind_f[ii][2], p);
        vector3 grad_lambda_j = grad_lambda(ind_f[ii][0]);
        vector3 grad_lambda_k = grad_lambda(ind_f[ii][1]);
        vector3 grad_lambda_l = grad_lambda(ind_f[ii][2]);
        return (lambda_l * grad_lambda_k + lambda_k * grad_lambda_l).cross(grad_lambda_j) -
               2.0 * (lambda_l * grad_lambda_j + lambda_j * grad_lambda_l).cross(grad_lambda_k) +
               (lambda_k * grad_lambda_j + lambda_j * grad_lambda_k).cross(grad_lambda_l);
    }
    // Второй полный
    else if(i < 30)
    {
        return vector3(0.0, 0.0, 0.0);
    }

    return vector3();
}

vector3 tetrahedron_base::kerw(size_t i, const point & p) const
{
    assert(i < basis->tet_ker_bf_num);
    if(i < 4)   return grad_lambda(i);
    if(i < 10)  return w(i + 2, p);
    if(i < 20)  return w(i + 10, p);
    return vector3();
}

void tetrahedron_base::init(const basis_type * basis)
{
    init_3d(basis->tet_int);
    this->basis = basis;
}

trio<double, vector3, cvector3>
tetrahedron_base::diff_normL2(const array_t<complex<double> > & q, eval_func func, void * data)
{
    const tetrahedron_integration * tet_integration = &(basis->tet_int);
    complex<double> result = 0.0;
    vector3 result_v3(0.0, 0.0, 0.0);
    cvector3 result_cv3(0.0, 0.0, 0.0);
    for(size_t k = 0; k < tet_integration->get_gauss_num(); k++)
    {
        cvector3 val(0.0, 0.0, 0.0);
        for(size_t i = 0; i < basis->tet_bf_num; i++)
            val = val + q[i] * cvector3(w(i, get_gauss_point(k)));
        cvector3 func_d = func(get_gauss_point(k), get_phys_area(), data) - val;
        result += tet_integration->get_gauss_weight(k) * func_d.norm2();
        for(size_t i = 0; i < 3; i++)
        {
            double re = func_d[i].real() * func_d[i].real();
            double im = func_d[i].imag() * func_d[i].imag();
            func_d[i] = complex<double>(re, im);
            result_v3[i] += tet_integration->get_gauss_weight(k) * (re + im);
        }
        result_cv3 += tet_integration->get_gauss_weight(k) * func_d;
    }
    result *= get_jacobian();
    result_v3 *= get_jacobian();
    result_cv3 *= get_jacobian();
    return make_trio(result.real(), result_v3, result_cv3);
}

trio<double, vector3, cvector3>
tetrahedron_base::diff_normL2(const array_t<complex<double> > & q, const array_t<complex<double> > & q_true)
{
    const tetrahedron_integration * tet_integration = &(basis->tet_int);
    complex<double> result = 0.0;
    vector3 result_v3(0.0, 0.0, 0.0);
    cvector3 result_cv3(0.0, 0.0, 0.0);
    for(size_t k = 0; k < tet_integration->get_gauss_num(); k++)
    {
        cvector3 val(0.0, 0.0, 0.0), val_true(0.0, 0.0, 0.0);
        for(size_t i = 0; i < basis->tet_bf_num; i++)
        {
            val = val + q[i] * cvector3(w(i, get_gauss_point(k)));
            val_true = val_true + q_true[i] * cvector3(w(i, get_gauss_point(k)));
        }
        cvector3 func_d = val_true - val;
        result += tet_integration->get_gauss_weight(k) * func_d.norm2();
        for(size_t i = 0; i < 3; i++)
        {
            double re = func_d[i].real() * func_d[i].real();
            double im = func_d[i].imag() * func_d[i].imag();
            func_d[i] = complex<double>(re, im);
            result_v3[i] += tet_integration->get_gauss_weight(k) * (re + im);
        }
        result_cv3 += tet_integration->get_gauss_weight(k) * func_d;
    }
    result *= get_jacobian();
    result_v3 *= get_jacobian();
    result_cv3 *= get_jacobian();
    return make_trio(result.real(), result_v3, result_cv3);
}

trio<double, vector3, cvector3>
tetrahedron_base::normL2(eval_func func, void * data)
{
    const tetrahedron_integration * tet_integration = &(basis->tet_int);
    complex<double> result = 0.0;
    vector3 result_v3(0.0, 0.0, 0.0);
    cvector3 result_cv3(0.0, 0.0, 0.0);
    for(size_t k = 0; k < tet_integration->get_gauss_num(); k++)
    {
        cvector3 func_d = func(get_gauss_point(k), get_phys_area(), data);
        result += tet_integration->get_gauss_weight(k) * func_d.norm2();
        for(size_t i = 0; i < 3; i++)
        {
            double re = func_d[i].real() * func_d[i].real();
            double im = func_d[i].imag() * func_d[i].imag();
            func_d[i] = complex<double>(re, im);
            result_v3[i] += tet_integration->get_gauss_weight(k) * (re + im);
        }
        result_cv3 += tet_integration->get_gauss_weight(k) * func_d;
    }
    result *= get_jacobian();
    result_v3 *= get_jacobian();
    result_cv3 *= get_jacobian();
    return make_trio(result.real(), result_v3, result_cv3);
}

trio<double, vector3, cvector3>
tetrahedron_base::normL2(const array_t<complex<double> > & q_true)
{
    const tetrahedron_integration * tet_integration = &(basis->tet_int);
    complex<double> result = 0.0;
    vector3 result_v3(0.0, 0.0, 0.0);
    cvector3 result_cv3(0.0, 0.0, 0.0);
    for(size_t k = 0; k < tet_integration->get_gauss_num(); k++)
    {
        cvector3 func_d;
        for(size_t i = 0; i < basis->tet_bf_num; i++)
            func_d = func_d + q_true[i] * cvector3(w(i, get_gauss_point(k)));
        result += tet_integration->get_gauss_weight(k) * func_d.norm2();
        for(size_t i = 0; i < 3; i++)
        {
            double re = func_d[i].real() * func_d[i].real();
            double im = func_d[i].imag() * func_d[i].imag();
            func_d[i] = complex<double>(re, im);
            result_v3[i] += tet_integration->get_gauss_weight(k) * (re + im);
        }
        result_cv3 += tet_integration->get_gauss_weight(k) * func_d;
    }
    result *= get_jacobian();
    result_v3 *= get_jacobian();
    result_cv3 *= get_jacobian();
    return make_trio(result.real(), result_v3, result_cv3);
}

// ============================================================================

// Локальная матрица полного пространства
matrix_t<complex<double> >
tetrahedron::MpG()
{
    const tetrahedron_integration * tet_integration = &(basis->tet_int);
    matrix_t<complex<double> > matr(basis->tet_bf_num, basis->tet_bf_num);
    phys_area * phys = get_phys_area_ptr();

    for(size_t i = 0; i < basis->tet_bf_num; i++)
        for(size_t j = 0; j < basis->tet_bf_num; j++)
            matr[i][j] = 0.0;

    for(size_t k = 0; k < tet_integration->get_gauss_num(); k++)
    {
        phys->sigma[threads_config::matrix_full].set_x(get_gauss_point(k).x);
        phys->sigma[threads_config::matrix_full].set_y(get_gauss_point(k).y);
        phys->sigma[threads_config::matrix_full].set_z(get_gauss_point(k).z);
        double sigma = 0.0;
        phys->sigma[threads_config::matrix_full].calculate(sigma);
        complex<double> k2(- phys->epsilon * phys->omega * phys->omega, phys->omega * sigma);

        for(size_t i = 0; i < basis->tet_bf_num; i++)
        {
            for(size_t j = 0; j <= i; j++)
            {
                // Интеграл от бф
                vector3 wi = w(i, get_gauss_point(k));
                vector3 wj = w(j, get_gauss_point(k));
                // Интеграл от ротора бф
                vector3 curlwi = rotw(i, get_gauss_point(k));
                vector3 curlwj = rotw(j, get_gauss_point(k));
                // Почти элемент локальной матрицы
                matr[i][j] += tet_integration->get_gauss_weight(k) * (wi * wj * k2 + curlwi * curlwj / phys->mu);
            }
        }
    }

    for(size_t i = 0; i < basis->tet_bf_num; i++)
    {
        for(size_t j = 0; j <= i; j++)
        {
            matr[i][j] *= get_jacobian();
            matr[j][i] = matr[i][j];
        }
    }
    return matr;
}

// Локальная правая часть
array_t<complex<double> >
tetrahedron::rp(eval_func func, void * data)
{
    const tetrahedron_integration * tet_integration = &(basis->tet_int);
    array_t<complex<double> > arr(basis->tet_bf_num);
    for(size_t i = 0; i < basis->tet_bf_num; i++)
    {
        complex<double> value(0, 0);
        for(size_t k = 0; k < tet_integration->get_gauss_num(); k++)
        {
            // Интеграл от ф-и правой части на бф
            vector3 wi = w(i, get_gauss_point(k));
            cvector3 f = func(get_gauss_point(k), get_phys_area(), data);
            // Почти элемент локальной правой части
            value += tet_integration->get_gauss_weight(k) * f * wi;
        }
        arr[i] = value * get_jacobian();
    }
    return arr;
}

// Локальная матрица ядра
matrix_t<complex<double> >
tetrahedron::K()
{
    const tetrahedron_integration * tet_integration = &(basis->tet_int);
    matrix_t<complex<double> > matr(basis->tet_ker_bf_num, basis->tet_ker_bf_num);
    phys_area * phys = get_phys_area_ptr();

    for(size_t i = 0; i < basis->tet_ker_bf_num; i++)
        for(size_t j = 0; j < basis->tet_ker_bf_num; j++)
            matr[i][j] = 0.0;

    for(size_t k = 0; k < tet_integration->get_gauss_num(); k++)
    {
        phys->sigma[threads_config::matrix_ker].set_x(get_gauss_point(k).x);
        phys->sigma[threads_config::matrix_ker].set_y(get_gauss_point(k).y);
        phys->sigma[threads_config::matrix_ker].set_z(get_gauss_point(k).z);
        double sigma = 0.0;
        phys->sigma[threads_config::matrix_ker].calculate(sigma);
        complex<double> k2(- phys->epsilon * phys->omega * phys->omega, phys->omega * sigma);

        for(size_t i = 0; i < basis->tet_ker_bf_num; i++)
        {
            for(size_t j = 0; j <= i; j++)
            {
                // Интегралы от базисных функций ядра
                vector3 kerwi = kerw(i, get_gauss_point(k));
                vector3 kerwj = kerw(j, get_gauss_point(k));
                // Почти элемент локальной матрицы
                matr[i][j] += tet_integration->get_gauss_weight(k) * kerwi * kerwj * k2;
            }
        }
    }

    for(size_t i = 0; i < basis->tet_ker_bf_num; i++)
    {
        for(size_t j = 0; j <= i; j++)
        {
            matr[i][j] *= get_jacobian();
            matr[j][i] = matr[i][j];
        }
    }
    return matr;
}

// ============================================================================
