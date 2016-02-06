#include "tetrahedron.h"

tetrahedron_pml::tetrahedron_pml()
{
    phys_pml = NULL;
    get_s = NULL;
    jacobian_pml = 0.0;
}

// Инициализация PML-координат
void tetrahedron_pml::init_pml(cvector3(* get_s)(const point &, const tetrahedron_pml *, const phys_pml_area *), const phys_pml_area * phys_pml, const cpoint * nodes_pml)
{
    const tet_integration_config * tet_integration = &(basis->tet_int);

    this->get_s = get_s;
    this->phys_pml = phys_pml;

    matrix_t<complex<double>, 4, 4> D;
    for(size_t i = 0; i < 3; i++)
        for(size_t j = 0; j < 4; j++)
            D[i][j] = nodes_pml[j][i];
    for(size_t j = 0; j < 4; j++)
        D[3][j] = 1.0;
    complex<double> D_det;
    L_pml = inverse(D, D_det);

    jacobian_pml = abs(D_det);

    // Перевод точек Гаусса с мастер-элемента на текущий тетраэдр
    gauss_points_pml.resize(tet_integration->gauss_num);
    for(size_t i = 0; i < 3; i++)
    {
        for(size_t j = 0 ; j < tet_integration->gauss_num; j++)
        {
            gauss_points_pml[j][i] = 0;
            for(size_t k = 0; k < 4; k++)
                gauss_points_pml[j][i] += D[i][k] * tet_integration->gauss_points_master[j][k];
        }
    }
}

// Получить указатель на обычный тетраэдр
const tetrahedron * tetrahedron_pml::to_std() const
{
    const tetrahedron_base * base = this;
    const tetrahedron * tet = static_cast<const tetrahedron *>(base);
    return tet;
}

// L-координаты (в PML)
complex<double> tetrahedron_pml::lambda_pml(size_t i, const cpoint & p) const
{
    assert(i < 4); // Если i = 4, то явно где-то косяк
    return L_pml[i][3] + L_pml[i][0] * p.x + L_pml[i][1] * p.y + L_pml[i][2] * p.z;
}

// Градиент L-координаты (в PML)
cvector3 tetrahedron_pml::grad_lambda_pml(size_t i) const
{
    assert(i < 4); // Если i = 4, то явно где-то косяк
    return cvector3(L_pml[i][0], L_pml[i][1], L_pml[i][2]);
}

// Базисные функции (в PML)
cvector3 tetrahedron_pml::w_pml(size_t i, const cpoint & p) const
{
    using namespace tet_basis_indexes;
    assert(i < basis->tet_bf_num);

    // Первый неполный
    if(i < 6)
    {
        return lambda_pml(ind_e[i][0], p) * grad_lambda_pml(ind_e[i][1]) -
               lambda_pml(ind_e[i][1], p) * grad_lambda_pml(ind_e[i][0]);
    }
    // Первый полный
    else if(i < 12)
    {
        size_t ii = i - 6;
        return lambda_pml(ind_e[ii][0], p) * grad_lambda_pml(ind_e[ii][1]) +
               lambda_pml(ind_e[ii][1], p) * grad_lambda_pml(ind_e[ii][0]);
    }
    // Второй неполный
    else if(i < 16)
    {
        size_t ii = i - 12;
        return lambda_pml(ind_f[ii][1], p) * lambda_pml(ind_f[ii][2], p) * grad_lambda_pml(ind_f[ii][0]) +
               lambda_pml(ind_f[ii][0], p) * lambda_pml(ind_f[ii][2], p) * grad_lambda_pml(ind_f[ii][1]) -
               2.0 * lambda_pml(ind_f[ii][0], p) * lambda_pml(ind_f[ii][1], p) * grad_lambda_pml(ind_f[ii][2]);
    }
    else if(i < 20)
    {
        size_t ii = i - 16;
        return lambda_pml(ind_f[ii][1], p) * lambda_pml(ind_f[ii][2], p) * grad_lambda_pml(ind_f[ii][0]) -
               2.0 * lambda_pml(ind_f[ii][0], p) * lambda_pml(ind_f[ii][2], p) * grad_lambda_pml(ind_f[ii][1]) +
               lambda_pml(ind_f[ii][0], p) * lambda_pml(ind_f[ii][1], p) * grad_lambda_pml(ind_f[ii][2]);
    }
    // Второй полный
    else if(i < 24)
    {
        size_t ii = i - 20;
        return lambda_pml(ind_f[ii][1], p) * lambda_pml(ind_f[ii][2], p) * grad_lambda_pml(ind_f[ii][0]) +
               lambda_pml(ind_f[ii][0], p) * lambda_pml(ind_f[ii][2], p) * grad_lambda_pml(ind_f[ii][1]) +
               lambda_pml(ind_f[ii][0], p) * lambda_pml(ind_f[ii][1], p) * grad_lambda_pml(ind_f[ii][2]);
    }
    else if(i < 30)
    {
        size_t ii = i - 24;
        return lambda_pml(ind_e[ii][1], p) * (2.0 * lambda_pml(ind_e[ii][0], p) - lambda_pml(ind_e[ii][1], p)) * grad_lambda_pml(ind_e[ii][0]) -
               lambda_pml(ind_e[ii][0], p) * (2.0 * lambda_pml(ind_e[ii][1], p) - lambda_pml(ind_e[ii][0], p)) * grad_lambda_pml(ind_e[ii][1]);
    }

    return cvector3();
}

// Роторы базисных функций (в PML)
cvector3 tetrahedron_pml::rotw_pml(size_t i, const cpoint & p, const point & p_non_PML) const
{
    using namespace tet_basis_indexes;
    assert(i < basis->tet_bf_num);
    assert(phys_pml);
    assert(get_s);

    // Первый неполный
    if(i < 6)
    {
        cvector3 grad1 = grad_lambda_pml(ind_e[i][0]);
        cvector3 grad2 = grad_lambda_pml(ind_e[i][1]);

        cvector3 s = get_s(p_non_PML, this, phys_pml);
        for(size_t ind = 0; ind < 3; ind++)
        {
            grad1[ind] /= s[ind];
            grad2[ind] /= s[ind];
        }
        return 2.0 * grad1.cross(grad2);
    }
    // Первый полный
    else if(i < 12)
    {
        return cvector3(0.0, 0.0, 0.0);
    }
    // Второй неполный
    else if(i < 16)
    {
        size_t ii = i - 12;
        complex<double> lambda_j = lambda_pml(ind_f[ii][0], p);
        complex<double> lambda_k = lambda_pml(ind_f[ii][1], p);
        complex<double> lambda_l = lambda_pml(ind_f[ii][2], p);
        cvector3 grad_lambda_j = grad_lambda_pml(ind_f[ii][0]);
        cvector3 grad_lambda_k = grad_lambda_pml(ind_f[ii][1]);
        cvector3 grad_lambda_l = grad_lambda_pml(ind_f[ii][2]);

        cvector3 s = get_s(p_non_PML, this, phys_pml);
        for(size_t ind = 0; ind < 3; ind++)
        {
            grad_lambda_j[ind] /= s[ind];
            grad_lambda_k[ind] /= s[ind];
            grad_lambda_l[ind] /= s[ind];
        }

        return (lambda_l * grad_lambda_k + lambda_k * grad_lambda_l).cross(grad_lambda_j) +
               (lambda_l * grad_lambda_j + lambda_j * grad_lambda_l).cross(grad_lambda_k) -
               2.0 * (lambda_k * grad_lambda_j + lambda_j * grad_lambda_k).cross(grad_lambda_l);
    }
    else if(i < 20)
    {
        size_t ii = i - 16;
        complex<double> lambda_j = lambda_pml(ind_f[ii][0], p);
        complex<double> lambda_k = lambda_pml(ind_f[ii][1], p);
        complex<double> lambda_l = lambda_pml(ind_f[ii][2], p);
        cvector3 grad_lambda_j = grad_lambda_pml(ind_f[ii][0]);
        cvector3 grad_lambda_k = grad_lambda_pml(ind_f[ii][1]);
        cvector3 grad_lambda_l = grad_lambda_pml(ind_f[ii][2]);

        cvector3 s = get_s(p_non_PML, this, phys_pml);
        for(size_t ind = 0; ind < 3; ind++)
        {
            grad_lambda_j[ind] /= s[ind];
            grad_lambda_k[ind] /= s[ind];
            grad_lambda_l[ind] /= s[ind];
        }

        return (lambda_l * grad_lambda_k + lambda_k * grad_lambda_l).cross(grad_lambda_j) -
               2.0 * (lambda_l * grad_lambda_j + lambda_j * grad_lambda_l).cross(grad_lambda_k) +
               (lambda_k * grad_lambda_j + lambda_j * grad_lambda_k).cross(grad_lambda_l);
    }
    // Второй полный
    else if(i < 30)
    {
        return cvector3(0.0, 0.0, 0.0);
    }

    return cvector3();
}

// Базисные функции ядра (в PML)
cvector3 tetrahedron_pml::kerw_pml(size_t i, const cpoint & p, const point & p_non_PML) const
{
    assert(i < basis->tet_ker_bf_num);
    assert(phys_pml);
    assert(get_s);
    if(i < 4)
    {
        cvector3 s = get_s(p_non_PML, this, phys_pml);
        cvector3 grad_lambda_i = grad_lambda(i);
        for(size_t ind = 0; ind < 3; ind++)
            grad_lambda_i[ind] /= s[ind];
        return grad_lambda_i;
    }
    if(i < 10)  return w_pml(i + 2, p);
    if(i < 20)  return w_pml(i + 10, p);
    return cvector3();
}

// Интеграл от бф
complex<double> tetrahedron_pml::integrate_w(size_t i, size_t j) const
{
//#if defined __GNUC__
//#warning conjugate
//#endif
    const tet_integration_config * tet_integration = &(basis->tet_int);
    complex<double> result = 0.0;
    for(size_t k = 0; k < tet_integration->gauss_num; k++)
        result += tet_integration->gauss_weights[k] *
                  w_pml(i, gauss_points_pml[k]).cj() *
                  w_pml(j, gauss_points_pml[k]).cj();
    return result * jacobian_pml;
}

// Интеграл от ротора бф
complex<double> tetrahedron_pml::integrate_rotw(size_t i, size_t j) const
{
//#if defined __GNUC__
//#warning conjugate
//#endif
    const tet_integration_config * tet_integration = &(basis->tet_int);
    complex<double> result = 0.0;
    for(size_t k = 0; k < tet_integration->gauss_num; k++)
        result += tet_integration->gauss_weights[k] *
                  rotw_pml(i, gauss_points_pml[k], gauss_points[k]).cj() *
                  rotw_pml(j, gauss_points_pml[k], gauss_points[k]).cj();
    return result * jacobian_pml;
}

// Интеграл от ф-и правой части на бф
complex<double> tetrahedron_pml::integrate_fw(eval_func func, size_t i, void * data) const
{
//#if defined __GNUC__
//#warning conjugate
//#endif
    const tet_integration_config * tet_integration = &(basis->tet_int);
    complex<double> result = 0.0;
    for(size_t k = 0; k < tet_integration->gauss_num; k++)
        result += tet_integration->gauss_weights[k] *
                  func(gauss_points[k], get_phys_area(), data) *
                  w_pml(i, gauss_points_pml[k]).cj();
    return result * jacobian_pml;
}

// Интегралы от базисных функций ядра
complex<double> tetrahedron_pml::integrate_kerw(size_t i, size_t j) const
{
//#if defined __GNUC__
//#warning conjugate
//#endif
    const tet_integration_config * tet_integration = &(basis->tet_int);
    complex<double> result = 0.0;
    for(size_t k = 0; k < tet_integration->gauss_num; k++)
        result += tet_integration->gauss_weights[k] *
                  kerw_pml(i, gauss_points_pml[k], gauss_points[k]).cj() *
                  kerw_pml(j, gauss_points_pml[k], gauss_points[k]).cj();
    return result * jacobian;
}

// Локальная матрица жескости
matrix_t<complex<double> >
tetrahedron_pml::G() const
{
    matrix_t<complex<double> > matr(basis->tet_bf_num, basis->tet_bf_num);
    for(size_t i = 0; i < basis->tet_bf_num; i++)
        for(size_t j = 0; j <= i; j++)
            matr[j][i] = matr[i][j] = integrate_rotw(i, j);
    return matr;
}

// Локальная матрица массы
matrix_t<complex<double> >
tetrahedron_pml::M() const
{
    matrix_t<complex<double> > matr(basis->tet_bf_num, basis->tet_bf_num);
    for(size_t i = 0; i < basis->tet_bf_num; i++)
        for(size_t j = 0; j <= i; j++)
            matr[j][i] = matr[i][j] = integrate_w(i, j);
    return matr;
}

// Локальная правая часть
array_t<complex<double> >
tetrahedron_pml::rp(eval_func func, void * data) const
{
    array_t<complex<double> > arr(basis->tet_bf_num);
    for(size_t i = 0; i < basis->tet_bf_num; i++)
        arr[i] = integrate_fw(func, i, data);
    return arr;
}

// Локальная матрица ядра
matrix_t<complex<double> >
tetrahedron_pml::K() const
{
    matrix_t<complex<double> > matr(basis->tet_ker_bf_num, basis->tet_ker_bf_num);
    for(size_t i = 0; i < basis->tet_ker_bf_num; i++)
        for(size_t j = 0; j <= i; j++)
            matr[j][i] = matr[i][j] = integrate_kerw(i, j);
    return matr;
}
