#include "tetrahedron.h"

// ============================================================================

tetrahedron_base::tetrahedron_base()
{
    for(size_t i = 0; i < 4; i++)
        nodes[i] = NULL;
    for(size_t i = 0; i < 6; i++)
        edges[i] = NULL;
#if BASIS_ORDER >= 2
    for(size_t i = 0; i < 4; i++)
        faces[i] = NULL;
#endif
    phys = NULL;
}

const point & tetrahedron_base::get_node(size_t i) const
{
    assert(i < 4);
    assert(nodes[i] != NULL);
    return (* nodes[i]);
}

const edge & tetrahedron_base::get_edge(size_t i) const
{
    assert(i < 6);
    assert(edges[i] != NULL);
    return (* edges[i]);
}

#if BASIS_ORDER >= 2
const face & tetrahedron_base::get_face(size_t i) const
{
    assert(i < 4);
    assert(faces[i] != NULL);
    return (* faces[i]);
}
#endif

const phys_area & tetrahedron_base::get_phys_area() const
{
    assert(phys != NULL);
    return * phys;
}

double tetrahedron_base::lambda(size_t i, const point & p) const
{
    assert(i < 4); // Если i == 4, то явно где-то косяк
    return L[i][3] + L[i][0] * p.x + L[i][1] * p.y + L[i][2] * p.z;
}

vector3 tetrahedron_base::grad_lambda(size_t i) const
{
    assert(i < 4); // Если i == 4, то явно где-то косяк
    return vector3(L[i][0], L[i][1], L[i][2]);
}

vector3 tetrahedron_base::w(size_t i, const point & p) const
{
    using namespace tet_basis_indexes;
    assert(i < basis::tet_bf_num);

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

vector3 tetrahedron_base::rotw(size_t i, const point & p) const
{
    using namespace tet_basis_indexes;
    assert(i < basis::tet_bf_num);
    MAYBE_UNUSED(p);

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

void tetrahedron_base::init()
{
    using namespace tet_integration;

    matrix_t<double, 4, 4> D;
    for(size_t i = 0; i < 3; i++)
        for(size_t j = 0; j < 4; j++)
            D[i][j] = get_node(j)[i];
    for(size_t j = 0; j < 4; j++)
        D[3][j] = 1.0;
    double D_det;
    L = inverse(D, D_det);

    jacobian = fabs(D_det);

    // Перевод точек Гаусса с мастер-элемента на текущий тетраэдр
    for(size_t i = 0; i < 3; i++)
    {
        for(size_t j = 0 ; j < gauss_num; j++)
        {
            gauss_points[j][i] = 0;
            for(size_t k = 0; k < 4; k++)
                gauss_points[j][i] += D[i][k] * gauss_points_master[j][k];
        }
    }

    // Расчет барицентра
    for(size_t i = 0; i < 3; i++)
    {
        barycenter[i] = 0.0;
        for(size_t j = 0; j < 4; j++)
            barycenter[i] += get_node(j)[i];
        barycenter[i] /= 4.0;
    }

    // Представление прямых тетраэдра в параметрическом виде a*t + b
    // Само ребо получается при 0<=t<=1
    for(size_t i = 0; i < 3; i++)
    {
        edges_a[0][i] = get_node(1)[i] - get_node(0)[i];
        edges_b[0][i] = get_node(0)[i];
        edges_a[1][i] = get_node(2)[i] - get_node(0)[i];
        edges_b[1][i] = get_node(0)[i];
        edges_a[2][i] = get_node(3)[i] - get_node(0)[i];
        edges_b[2][i] = get_node(0)[i];
        edges_a[3][i] = get_node(2)[i] - get_node(1)[i];
        edges_b[3][i] = get_node(1)[i];
        edges_a[4][i] = get_node(3)[i] - get_node(1)[i];
        edges_b[4][i] = get_node(1)[i];
        edges_a[5][i] = get_node(3)[i] - get_node(2)[i];
        edges_b[5][i] = get_node(2)[i];
    }
}

bool tetrahedron_base::inside(const point & p) const
{
    for(size_t i = 0; i < 4; i++)
    {
        double a = lambda(i, p);
        if(a -1e-12 > 1.0 || a +1e-12 < 0.0)
            return false;
    }
    return true;
}

bool tetrahedron_base::inside(double x, double y, double z) const
{
    return inside(point(x, y, z));
}

bool tetrahedron_base::inside_tree(double x0, double x1, double y0, double y1, double z0, double z1) const
{
    // Тетраэдр внутри куба
    if(barycenter.inside(x0, x1, y0, y1, z0, z1))
        return true;
    for(size_t i = 0; i < 4; i++)
        if(get_node(i).inside(x0, x1, y0, y1, z0, z1))
            return true;

    // Куб внутри тетраэдра
    if(inside(x0, y0, z0) || inside(x1, y0, z0) || inside(x0, y1, z0) || inside(x1, y1, z0) ||
            inside(x0, y0, z1) || inside(x1, y0, z1) || inside(x0, y1, z1) || inside(x1, y1, z1))
        return true;

    // Пересечение куба и тетраэдра, не обработанное выше
    double t[6][6];     // Параметры прямых для всех прямых (6) и всех плоскостей (6)
    double cords_t[6][6][3];    // Прямая, плоскость, координата
    for(size_t i = 0; i < 6; i++)
    {
        t[i][0] = (x0 - edges_b[i][0]) / edges_a[i][0]; // x = x0
        t[i][1] = (x1 - edges_b[i][0]) / edges_a[i][0]; // x = x1
        t[i][2] = (y0 - edges_b[i][1]) / edges_a[i][1]; // y = y0
        t[i][3] = (y1 - edges_b[i][1]) / edges_a[i][1]; // y = y1
        t[i][4] = (z0 - edges_b[i][2]) / edges_a[i][2]; // z = z0
        t[i][5] = (z1 - edges_b[i][2]) / edges_a[i][2]; // z = z1
    }

    for(size_t i = 0; i < 6; i++)
        for(size_t j = 0; j < 6; j++)
            for(size_t k = 0; k < 3; k++)
                cords_t[i][j][k] = edges_a[i][k] * t[i][j] + edges_b[i][k];

    for(size_t i = 0; i < 6; i++)   // Берем прямую и проверяем, что пересечение с плоскостью попадает в раcсматриваемый отрезок прямой и в рассатриваемую часть плоскости
    {
        if(
            (t[i][0] >= 0.0 && t[i][0] <= 1.0 && cords_t[i][0][1] >= y0 && cords_t[i][0][1] >= y1 && cords_t[i][0][2] >= z0 && cords_t[i][0][2] <= z1) || // x = x0
            (t[i][1] >= 0.0 && t[i][1] <= 1.0 && cords_t[i][1][1] >= y0 && cords_t[i][1][1] >= y1 && cords_t[i][1][2] >= z0 && cords_t[i][1][2] <= z1) || // x = x1
            (t[i][2] >= 0.0 && t[i][2] <= 1.0 && cords_t[i][2][0] >= x0 && cords_t[i][2][0] <= x1 && cords_t[i][2][2] >= z0 && cords_t[i][2][2] <= z1) || // y = y0
            (t[i][3] >= 0.0 && t[i][3] <= 1.0 && cords_t[i][3][0] >= x0 && cords_t[i][3][0] <= x1 && cords_t[i][3][2] >= z0 && cords_t[i][3][2] <= z1) || // y = y1
            (t[i][4] >= 0.0 && t[i][4] <= 1.0 && cords_t[i][4][0] >= x0 && cords_t[i][4][0] <= x1 && cords_t[i][4][1] >= y0 && cords_t[i][4][1] <= y1) || // z = z0
            (t[i][5] >= 0.0 && t[i][5] <= 1.0 && cords_t[i][5][0] >= x0 && cords_t[i][5][0] <= x1 && cords_t[i][5][1] >= y0 && cords_t[i][5][1] <= y1)    // z = z1
        )
            return true;
    }
    return false;
}

double tetrahedron_base::diff_normL2(const array_t<complex<double>, basis::tet_bf_num> & q, cvector3(*func)(const point &)) const
{
    using namespace tet_integration;
    using namespace basis;
    complex<double> result = 0.0;
    for(size_t k = 0; k < gauss_num; k++)
    {
        cvector3 val(0.0, 0.0, 0.0);
        for(size_t i = 0; i < tet_bf_num; i++)
            val = val + q[i] * cvector3(w(i, gauss_points[k]));
        cvector3 func_d = func(gauss_points[k]) - val;
        //result += gauss_weights[k] * (func_d * func_d.cj());
        result += gauss_weights[k] * func_d.norm2();
    }
    result *= jacobian;
    return result.real();
}

double tetrahedron_base::diff_normL2(const array_t<complex<double>, basis::tet_bf_num> & q, const array_t<complex<double>, basis::tet_bf_num> & q_true) const
{
    using namespace tet_integration;
    using namespace basis;
    complex<double> result = 0.0;
    for(size_t k = 0; k < gauss_num; k++)
    {
        cvector3 val(0.0, 0.0, 0.0), val_true(0.0, 0.0, 0.0);
        for(size_t i = 0; i < tet_bf_num; i++)
        {
            val = val + q[i] * cvector3(w(i, gauss_points[k]));
            val_true = val_true + q_true[i] * cvector3(w(i, gauss_points[k]));
        }
        cvector3 func_d = val_true - val;
        //result += gauss_weights[k] * (func_d * func_d.cj());
        result += gauss_weights[k] * func_d.norm2();
    }
    result *= jacobian;
    return result.real();
}

double tetrahedron_base::normL2(cvector3(*func)(const point &)) const
{
    using namespace tet_integration;
    complex<double> result = 0.0;
    for(size_t k = 0; k < gauss_num; k++)
    {
        cvector3 func_d = func(gauss_points[k]);
        //result += gauss_weights[k] * (func_d * func_d.cj());
        result += gauss_weights[k] * func_d.norm2();
    }
    result *= jacobian;
    return result.real();
}

double tetrahedron_base::normL2(const array_t<complex<double>, basis::tet_bf_num> & q_true) const
{
    using namespace tet_integration;
    using namespace basis;
    complex<double> result = 0.0;
    for(size_t k = 0; k < gauss_num; k++)
    {
        cvector3 func_d;
        for(size_t i = 0; i < tet_bf_num; i++)
            func_d = func_d + q_true[i] * cvector3(w(i, gauss_points[k]));
        //result += gauss_weights[k] * (func_d * func_d.cj());
        result += gauss_weights[k] * func_d.norm2();
    }
    result *= jacobian;
    return result.real();
}

// Локальная матрица ядра
matrix_t<double, 10, 10> tetrahedron_base::K() const
{
    matrix_t<double, 10, 10> matr;
    for(size_t i = 0; i < 10; i++)
        for(size_t j = 0; j <= i; j++)
            matr[j][i] = matr[i][j] = integrate_kerw(i, j);
    return matr;
}

#if defined __GNUC__
#warning "TODO"
#endif
// Локальная матрица проектора
matrix_t<double, 10, 12> tetrahedron_base::R() const
{
    // P - матрица инцидентности ребер и узлов, элемент Pij равен -1, если из i-го
    // узла выходит j-е ребро, 1, если ребро входит и 0 в остальных случаях
    static const double P[4][6] =
    {
        {-1.0, -1.0, -1.0,  0.0,  0.0,  0.0},
        { 1.0,  0.0,  0.0, -1.0, -1.0,  0.0},
        { 0.0,  1.0,  0.0,  1.0,  0.0, -1.0},
        { 0.0,  0.0,  1.0,  0.0,  1.0,  1.0}
    };
    // Матрица R имеет вид: | P  0 |
    //                      | 0  I |
    matrix_t<double, 10, 12> matr;
    for(size_t i = 0; i < 10; i++)
        for(size_t j = 0; j < 12; j++)
            matr[i][j] = 0.0;
    for(size_t i = 0; i < 4; i++)
        for(size_t j = 0; j < 6; j++)
            matr[i][j] = P[i][j];
    for(int i = 0; i < 6; i++)
        matr[i + 4][i + 6] = 1.0;
    return matr;
}

#if defined __GNUC__
#warning "TODO"
#endif
// Базисные функции ядра
vector3 tetrahedron_base::kerw(size_t i, const point & p) const
{
    if(i < 4)
        return grad_lambda(i);
    if(i < 10)
        return w(i + 2, p);
    cerr << "Error: Incorrect ker basis function number!" << endl;
    return vector3();
}

// Интегралы от базисных функций ядра
double tetrahedron_base::integrate_kerw(size_t i, size_t j) const
{
    using namespace tet_integration;
    double result = 0.0;
    for(size_t k = 0; k < gauss_num; k++)
        result += gauss_weights[k] *
                  kerw(i, gauss_points[k]) *
                  kerw(j, gauss_points[k]);
    return result * jacobian;
}

// ============================================================================

double tetrahedron::integrate_w(size_t i, size_t j) const
{
    using namespace tet_integration;
    double result = 0.0;
    for(size_t k = 0; k < gauss_num; k++)
        result += gauss_weights[k] *
                  w(i, gauss_points[k]) *
                  w(j, gauss_points[k]);
    return result * jacobian;
}

double tetrahedron::integrate_rotw(size_t i, size_t j) const
{
    using namespace tet_integration;
    double result = 0.0;
    for(size_t k = 0; k < gauss_num; k++)
        result += gauss_weights[k] *
                  rotw(i, gauss_points[k]) *
                  rotw(j, gauss_points[k]);
    return result * jacobian;
}

complex<double> tetrahedron::integrate_fw(cvector3(*func)(const point &, const phys_area &), size_t i) const
{
    using namespace tet_integration;
    complex<double> result = 0.0;
    for(size_t k = 0; k < gauss_num; k++)
        result += gauss_weights[k] *
                  func(gauss_points[k], get_phys_area()) *
                  w(i, gauss_points[k]);
    return result * jacobian;
}

matrix_t<double, basis::tet_bf_num, basis::tet_bf_num>
tetrahedron::G() const
{
    using namespace basis;
    matrix_t<double, tet_bf_num, tet_bf_num> matr;
    for(size_t i = 0; i < tet_bf_num; i++)
        for(size_t j = 0; j <= i; j++)
            matr[j][i] = matr[i][j] = integrate_rotw(i, j);
    return matr;
}

matrix_t<double, basis::tet_bf_num, basis::tet_bf_num>
tetrahedron::M() const
{
    using namespace basis;
    matrix_t<double, tet_bf_num, tet_bf_num> matr;
    for(size_t i = 0; i < tet_bf_num; i++)
        for(size_t j = 0; j <= i; j++)
            matr[j][i] = matr[i][j] = integrate_w(i, j);
    return matr;
}

array_t<complex<double>, basis::tet_bf_num>
tetrahedron::rp(cvector3(*func)(const point &, const phys_area &)) const
{
    using namespace basis;
    array_t<complex<double>, tet_bf_num> arr;
    for(size_t i = 0; i < tet_bf_num; i++)
        arr[i] = integrate_fw(func, i);
    return arr;
}

// ============================================================================
