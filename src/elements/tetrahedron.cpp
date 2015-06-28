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
matrix_t<double, 10, 12> tetrahedron_base::P() const
{
    // P11 - матрица инцидентности ребер и узлов, элемент P11ij равен -1, если из i-го
    // узла выходит j-е ребро, 1, если ребро входит и 0 в остальных случаях
    static const double P11[4][6] =
    {
        {-1.0, -1.0, -1.0,  0.0,  0.0,  0.0},
        { 1.0,  0.0,  0.0, -1.0, -1.0,  0.0},
        { 0.0,  1.0,  0.0,  1.0,  0.0, -1.0},
        { 0.0,  0.0,  1.0,  0.0,  1.0,  1.0}
    };
    // Матрица P имеет вид: | P11  0  |
    //                      | 0    I6 |
    matrix_t<double, 10, 12> matr;
    for(size_t i = 0; i < 10; i++)
        for(size_t j = 0; j < 12; j++)
            matr[i][j] = 0.0;
    for(size_t i = 0; i < 4; i++)
        for(size_t j = 0; j < 6; j++)
            matr[i][j] = P11[i][j];
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

matrix_t<double, 10, 12> tetrahedron_base::GetNodalStiff()
{
    double x[4], y[4], z[4];
    x[0] = nodes[0]->x;
    y[0] = nodes[0]->y;
    z[0] = nodes[0]->z;

    x[1] = nodes[1]->x;
    y[1] = nodes[1]->y;
    z[1] = nodes[1]->z;

    x[2] = nodes[2]->x;
    y[2] = nodes[2]->y;
    z[2] = nodes[2]->z;

    x[3] = nodes[3]->x;
    y[3] = nodes[3]->y;
    z[3] = nodes[3]->z;

    //определим коэффициенты
    double a1 = y[1]*z[2] - y[1]*z[3] - z[1]*y[2] + z[1]*y[3] + y[2]*z[3] - y[3]*z[2];
    double a2 = -y[0]*z[2]+y[0]*z[3]+z[0]*y[2]-z[0]*y[3]-y[2]*z[3]+y[3]*z[2];
    double a3 = y[0]*z[1]-y[0]*z[3]-z[0]*y[1]+z[0]*y[3]+y[1]*z[3]-z[1]*y[3];
    double a4 = -y[0]*z[1]+y[0]*z[2]+z[0]*y[1]-z[0]*y[2]-y[1]*z[2]+z[1]*y[2];

    double b1 = -x[1]*z[2]+x[1]*z[3]+z[1]*x[2]-z[1]*x[3]-x[2]*z[3]+x[3]*z[2];
    double b2 = x[0]*z[2]-x[0]*z[3]-z[0]*x[2]+z[0]*x[3]+x[2]*z[3]-x[3]*z[2];
    double b3 = -x[0]*z[1]+x[0]*z[3]+z[0]*x[1]-z[0]*x[3]-x[1]*z[3]+z[1]*x[3];
    double b4 = x[0]*z[1]-x[0]*z[2]-z[0]*x[1]+z[0]*x[2]+x[1]*z[2]-z[1]*x[2];

    double c1 = x[1]*y[2]-x[1]*y[3]-y[1]*x[2]+y[1]*x[3]+x[2]*y[3]-x[3]*y[2];
    double c2 = -x[0]*y[2]+x[0]*y[3]+y[0]*x[2]-y[0]*x[3]-x[2]*y[3]+x[3]*y[2];
    double c3 = x[0]*y[1]-x[0]*y[3]-y[0]*x[1]+y[0]*x[3]+x[1]*y[3]-y[1]*x[3];
    double c4 = -x[0]*y[1]+x[0]*y[2]+y[0]*x[1]-y[0]*x[2]-x[1]*y[2]+y[1]*x[2];

    double rdet =  x[0]*y[1]*z[2]-x[0]*y[1]*z[3]-x[0]*z[1]*y[2]+x[0]*z[1]*y[3]+x[0]*y
        [2]*z[3]-x[0]*y[3]*z[2]-y[0]*x[1]*z[2]+y[0]*x[1]*z[3]+y[0]*z[1]*x[2]-y[0]*z[1]*
        x[3]-y[0]*x[2]*z[3]+y[0]*x[3]*z[2]+z[0]*x[1]*y[2]-z[0]*x[1]*y[3]-z[0]*y[1]*x[2]
    +z[0]*y[1]*x[3]+z[0]*x[2]*y[3]-z[0]*x[3]*y[2]-x[1]*y[2]*z[3]+x[1]*y[3]*z[2]+y
        [1]*x[2]*z[3]-y[1]*x[3]*z[2]-z[1]*x[2]*y[3]+z[1]*x[3]*y[2];

    //определитель берем положительный
    rdet = jacobian;

    double C_G[20][20];

    C_G[0][0] = (a1*a1+b1*b1+c1*c1)/rdet/6.0;
    C_G[0][1] = (a2*a1+b2*b1+c2*c1)/rdet/6.0;
    C_G[0][2] = (a3*a1+b3*b1+c3*c1)/rdet/6.0;
    C_G[0][3] = (a1*a4+b1*b4+c1*c4)/rdet/6.0;
    C_G[0][4] = (a1*a1+a2*a1+b1*b1+b2*b1+c1*c1+c2*c1)/rdet/24.0;
    C_G[0][5] = (a1*a1+a3*a1+b1*b1+b3*b1+c1*c1+c3*c1)/rdet/24.0;
    C_G[0][6] = (a1*a1+a1*a4+b1*b1+b1*b4+c1*c1+c1*c4)/rdet/24.0;
    C_G[0][7] = (a2*a1+a3*a1+b2*b1+b3*b1+c2*c1+c3*c1)/rdet/24.0;
    C_G[0][8] = (a2*a1+a1*a4+b2*b1+b1*b4+c2*c1+c1*c4)/rdet/24.0;
    C_G[0][9] = (a3*a1+a1*a4+b3*b1+b1*b4+c3*c1+c1*c4)/rdet/24.0;
    C_G[0][10] = (a2*a1+a3*a1+a1*a1+b2*b1+b3*b1+b1*b1+c2*c1+c3*c1+c1*c1)/rdet
        /120.0;
    C_G[0][11] = (a2*a1+a1*a4+a1*a1+b2*b1+b1*b4+b1*b1+c2*c1+c1*c4+c1*c1)/rdet
        /120.0;
    C_G[0][12] = (a3*a1+a1*a4+a1*a1+b3*b1+b1*b4+b1*b1+c3*c1+c1*c4+c1*c1)/rdet
        /120.0;
    C_G[0][13] = (a3*a1+a1*a4+a2*a1+b3*b1+b1*b4+b2*b1+c3*c1+c1*c4+c2*c1)/rdet
        /120.0;
    C_G[0][14] = 0.0;
    C_G[0][15] = 0.0;
    C_G[0][16] = 0.0;
    C_G[0][17] = 0.0;
    C_G[0][18] = 0.0;
    C_G[0][19] = 0.0;
    C_G[1][0] = (a2*a1+b2*b1+c2*c1)/rdet/6.0;
    C_G[1][1] = (a2*a2+b2*b2+c2*c2)/rdet/6.0;
    C_G[1][2] = (b3*b2+a3*a2+c3*c2)/rdet/6.0;
    C_G[1][3] = (a2*a4+b2*b4+c2*c4)/rdet/6.0;
    C_G[1][4] = (a2*a1+a2*a2+b2*b1+b2*b2+c2*c1+c2*c2)/rdet/24.0;
    C_G[1][5] = (a2*a1+a3*a2+b2*b1+b3*b2+c2*c1+c3*c2)/rdet/24.0;
    C_G[1][6] = (a2*a1+a2*a4+b2*b1+b2*b4+c2*c1+c2*c4)/rdet/24.0;
    C_G[1][7] = (a2*a2+a3*a2+b2*b2+b3*b2+c2*c2+c3*c2)/rdet/24.0;
    C_G[1][8] = (a2*a2+a2*a4+b2*b2+b2*b4+c2*c2+c2*c4)/rdet/24.0;
    C_G[1][9] = (a3*a2+a2*a4+b3*b2+b2*b4+c3*c2+c2*c4)/rdet/24.0;
    C_G[1][10] = (a2*a2+a3*a2+a2*a1+b2*b2+b3*b2+b2*b1+c2*c2+c3*c2+c2*c1)/rdet
        /120.0;
    C_G[1][11] = (a2*a2+a2*a4+a2*a1+b2*b2+b2*b4+b2*b1+c2*c2+c2*c4+c2*c1)/rdet
        /120.0;
    C_G[1][12] = (a3*a2+a2*a4+a2*a1+b3*b2+b2*b4+b2*b1+c3*c2+c2*c4+c2*c1)/rdet
        /120.0;
    C_G[1][13] = (a3*a2+a2*a4+a2*a2+b3*b2+b2*b4+b2*b2+c3*c2+c2*c4+c2*c2)/rdet
        /120.0;
    C_G[1][14] = 0.0;
    C_G[1][15] = 0.0;
    C_G[1][16] = 0.0;
    C_G[1][17] = 0.0;
    C_G[1][18] = 0.0;
    C_G[1][19] = 0.0;
    C_G[2][0] = (a3*a1+b3*b1+c3*c1)/rdet/6.0;
    C_G[2][1] = (b3*b2+a3*a2+c3*c2)/rdet/6.0;
    C_G[2][2] = (a3*a3+b3*b3+c3*c3)/rdet/6.0;
    C_G[2][3] = (a3*a4+b3*b4+c3*c4)/rdet/6.0;
    C_G[2][4] = (a3*a1+a3*a2+b3*b1+b3*b2+c3*c1+c3*c2)/rdet/24.0;
    C_G[2][5] = (a3*a1+a3*a3+b3*b1+b3*b3+c3*c1+c3*c3)/rdet/24.0;
    C_G[2][6] = (a3*a1+a3*a4+b3*b1+b3*b4+c3*c1+c3*c4)/rdet/24.0;
    C_G[2][7] = (a3*a2+a3*a3+b3*b2+b3*b3+c3*c2+c3*c3)/rdet/24.0;
    C_G[2][8] = (a3*a2+a3*a4+b3*b2+b3*b4+c3*c2+c3*c4)/rdet/24.0;
    C_G[2][9] = (a3*a3+a3*a4+b3*b3+b3*b4+c3*c3+c3*c4)/rdet/24.0;
    C_G[2][10] = (a3*a2+a3*a3+a3*a1+b3*b2+b3*b3+b3*b1+c3*c2+c3*c3+c3*c1)/rdet
        /120.0;
    C_G[2][11] = (a3*a2+a3*a4+a3*a1+b3*b2+b3*b4+b3*b1+c3*c2+c3*c4+c3*c1)/rdet
        /120.0;
    C_G[2][12] = (a3*a3+a3*a4+a3*a1+b3*b3+b3*b4+b3*b1+c3*c3+c3*c4+c3*c1)/rdet
        /120.0;
    C_G[2][13] = (a3*a3+a3*a4+a3*a2+b3*b3+b3*b4+b3*b2+c3*c3+c3*c4+c3*c2)/rdet
        /120.0;
    C_G[2][14] = 0.0;
    C_G[2][15] = 0.0;
    C_G[2][16] = 0.0;
    C_G[2][17] = 0.0;
    C_G[2][18] = 0.0;
    C_G[2][19] = 0.0;
    C_G[3][0] = (a1*a4+b1*b4+c1*c4)/rdet/6.0;
    C_G[3][1] = (a2*a4+b2*b4+c2*c4)/rdet/6.0;
    C_G[3][2] = (a3*a4+b3*b4+c3*c4)/rdet/6.0;
    C_G[3][3] = (a4*a4+b4*b4+c4*c4)/rdet/6.0;
    C_G[3][4] = (a1*a4+a2*a4+b1*b4+b2*b4+c1*c4+c2*c4)/rdet/24.0;
    C_G[3][5] = (a3*a4+a1*a4+b3*b4+b1*b4+c3*c4+c1*c4)/rdet/24.0;
    C_G[3][6] = (a4*a4+a1*a4+b4*b4+b1*b4+c4*c4+c1*c4)/rdet/24.0;
    C_G[3][7] = (a2*a4+a3*a4+b2*b4+b3*b4+c2*c4+c3*c4)/rdet/24.0;
    C_G[3][8] = (a2*a4+a4*a4+b2*b4+b4*b4+c2*c4+c4*c4)/rdet/24.0;
    C_G[3][9] = (a3*a4+a4*a4+b3*b4+b4*b4+c3*c4+c4*c4)/rdet/24.0;
    C_G[3][10] = (a2*a4+a3*a4+a1*a4+b2*b4+b3*b4+b1*b4+c2*c4+c3*c4+c1*c4)/rdet
        /120.0;
    C_G[3][11] = (a2*a4+a4*a4+a1*a4+b2*b4+b4*b4+b1*b4+c2*c4+c4*c4+c1*c4)/rdet
        /120.0;
    C_G[3][12] = (a3*a4+a4*a4+a1*a4+b3*b4+b4*b4+b1*b4+c3*c4+c4*c4+c1*c4)/rdet
        /120.0;
    C_G[3][13] = (a3*a4+a4*a4+a2*a4+b3*b4+b4*b4+b2*b4+c3*c4+c4*c4+c2*c4)/rdet
        /120.0;
    C_G[3][14] = 0.0;
    C_G[3][15] = 0.0;
    C_G[3][16] = 0.0;
    C_G[3][17] = 0.0;
    C_G[3][18] = 0.0;
    C_G[3][19] = 0.0;
    C_G[4][0] = (a1*a1+a2*a1+b1*b1+b2*b1+c1*c1+c2*c1)/rdet/24.0;
    C_G[4][1] = (a2*a1+a2*a2+b2*b1+b2*b2+c2*c1+c2*c2)/rdet/24.0;
    C_G[4][2] = (a3*a1+a3*a2+b3*b1+b3*b2+c3*c1+c3*c2)/rdet/24.0;
    C_G[4][3] = (a1*a4+a2*a4+b1*b4+b2*b4+c1*c4+c2*c4)/rdet/24.0;
    C_G[4][4] = (a1*a1+a2*a1+a2*a2+b1*b1+b2*b1+b2*b2+c1*c1+c2*c1+c2*c2)/rdet/
        60.0;
    C_G[4][5] = (a3*a1+a1*a1+2.0*a3*a2+a2*a1+b3*b1+b1*b1+2.0*b3*b2+b2*b1+c3*
        c1+c1*c1+2.0*c3*c2+c2*c1)/rdet/120.0;
    C_G[4][6] = (a1*a4+a1*a1+2.0*a2*a4+a2*a1+b1*b4+b1*b1+2.0*b2*b4+b2*b1+c1*
        c4+c1*c1+2.0*c2*c4+c2*c1)/rdet/120.0;
    C_G[4][7] = (a2*a1+2.0*a3*a1+a2*a2+a3*a2+b2*b1+2.0*b3*b1+b2*b2+b3*b2+c2*
        c1+2.0*c3*c1+c2*c2+c3*c2)/rdet/120.0;
    C_G[4][8] = (a2*a1+2.0*a1*a4+a2*a2+a2*a4+b2*b1+2.0*b1*b4+b2*b2+b2*b4+c2*
        c1+2.0*c1*c4+c2*c2+c2*c4)/rdet/120.0;
    C_G[4][9] = (a3*a1+a1*a4+a3*a2+a2*a4+b3*b1+b1*b4+b3*b2+b2*b4+c3*c1+c1*c4+
        c3*c2+c2*c4)/rdet/120.0;
    C_G[4][10] = (a2*a1+a3*a1+a1*a1+a2*a2+a3*a2+b2*b1+b3*b1+b1*b1+b2*b2+b3*b2
        +c2*c1+c3*c1+c1*c1+c2*c2+c3*c2)/rdet/360.0;
    C_G[4][11] = (a2*a1+a1*a4+a1*a1+a2*a2+a2*a4+b2*b1+b1*b4+b1*b1+b2*b2+b2*b4
        +c2*c1+c1*c4+c1*c1+c2*c2+c2*c4)/rdet/360.0;
    C_G[4][12] = (a3*a1+a1*a4+a1*a1+2.0*a3*a2+2.0*a2*a4+a2*a1+b3*b1+b1*b4+b1*
        b1+2.0*b3*b2+2.0*b2*b4+b2*b1+c3*c1+c1*c4+c1*c1+2.0*c3*c2+2.0*c2*c4+c2*c1)/rdet/
        720.0;
    C_G[4][13] = (2.0*a3*a1+2.0*a1*a4+a2*a1+a3*a2+a2*a4+a2*a2+2.0*b3*b1+2.0*
        b1*b4+b2*b1+b3*b2+b2*b4+b2*b2+2.0*c3*c1+2.0*c1*c4+c2*c1+c3*c2+c2*c4+c2*c2)/rdet
        /720.0;
    C_G[4][14] = -(a1*a1-a2*a2+b1*b1-b2*b2+c1*c1-c2*c2)/rdet/360.0;
    C_G[4][15] = (a2*a1+a3*a2+b2*b1+b3*b2+c2*c1+c3*c2)/rdet/360.0;
    C_G[4][16] = (a2*a1+a2*a4+b2*b1+b2*b4+c2*c1+c2*c4)/rdet/360.0;
    C_G[4][17] = (a2*a1+a3*a1+b2*b1+b3*b1+c2*c1+c3*c1)/rdet/360.0;
    C_G[4][18] = (a2*a1+a1*a4+b2*b1+b1*b4+c2*c1+c1*c4)/rdet/360.0;
    C_G[4][19] = 0.0;
    C_G[5][0] = (a1*a1+a3*a1+b1*b1+b3*b1+c1*c1+c3*c1)/rdet/24.0;
    C_G[5][1] = (a2*a1+a3*a2+b2*b1+b3*b2+c2*c1+c3*c2)/rdet/24.0;
    C_G[5][2] = (a3*a1+a3*a3+b3*b1+b3*b3+c3*c1+c3*c3)/rdet/24.0;
    C_G[5][3] = (a3*a4+a1*a4+b3*b4+b1*b4+c3*c4+c1*c4)/rdet/24.0;
    C_G[5][4] = (a3*a1+a1*a1+2.0*a3*a2+a2*a1+b3*b1+b1*b1+2.0*b3*b2+b2*b1+c3*
        c1+c1*c1+2.0*c3*c2+c2*c1)/rdet/120.0;
    C_G[5][5] = (a3*a3+a3*a1+a1*a1+b3*b3+b3*b1+b1*b1+c3*c3+c3*c1+c1*c1)/rdet/
        60.0;
    C_G[5][6] = (2.0*a3*a4+a3*a1+a1*a4+a1*a1+2.0*b3*b4+b3*b1+b1*b4+b1*b1+2.0*
        c3*c4+c3*c1+c1*c4+c1*c1)/rdet/120.0;
    C_G[5][7] = (a3*a2+a3*a3+2.0*a2*a1+a3*a1+b3*b2+b3*b3+2.0*b2*b1+b3*b1+c3*
        c2+c3*c3+2.0*c2*c1+c3*c1)/rdet/120.0;
    C_G[5][8] = (a3*a2+a3*a4+a2*a1+a1*a4+b3*b2+b3*b4+b2*b1+b1*b4+c3*c2+c3*c4+
        c2*c1+c1*c4)/rdet/120.0;
    C_G[5][9] = (a3*a3+a3*a4+a3*a1+2.0*a1*a4+b3*b3+b3*b4+b3*b1+2.0*b1*b4+c3*
        c3+c3*c4+c3*c1+2.0*c1*c4)/rdet/120.0;
    C_G[5][10] = (a3*a2+a3*a3+a3*a1+a2*a1+a1*a1+b3*b2+b3*b3+b3*b1+b2*b1+b1*b1
        +c3*c2+c3*c3+c3*c1+c2*c1+c1*c1)/rdet/360.0;
    C_G[5][11] = (2.0*a3*a2+2.0*a3*a4+a3*a1+a2*a1+a1*a4+a1*a1+2.0*b3*b2+2.0*
        b3*b4+b3*b1+b2*b1+b1*b4+b1*b1+2.0*c3*c2+2.0*c3*c4+c3*c1+c2*c1+c1*c4+c1*c1)/rdet
        /720.0;
    C_G[5][12] = (a3*a3+a3*a4+a3*a1+a1*a4+a1*a1+b3*b3+b3*b4+b3*b1+b1*b4+b1*b1
        +c3*c3+c3*c4+c3*c1+c1*c4+c1*c1)/rdet/360.0;
    C_G[5][13] = (a3*a3+a3*a4+a3*a2+a3*a1+2.0*a1*a4+2.0*a2*a1+b3*b3+b3*b4+b3*
        b2+b3*b1+2.0*b1*b4+2.0*b2*b1+c3*c3+c3*c4+c3*c2+c3*c1+2.0*c1*c4+2.0*c2*c1)/rdet/
        720.0;
    C_G[5][14] = (a3*a1+a3*a2+b3*b1+b3*b2+c3*c1+c3*c2)/rdet/360.0;
    C_G[5][15] = -(-a3*a3+a1*a1-b3*b3+b1*b1-c3*c3+c1*c1)/rdet/360.0;
    C_G[5][16] = (a3*a1+a3*a4+b3*b1+b3*b4+c3*c1+c3*c4)/rdet/360.0;
    C_G[5][17] = -(a2*a1+a3*a1+b2*b1+b3*b1+c2*c1+c3*c1)/rdet/360.0;
    C_G[5][18] = 0.0;
    C_G[5][19] = (a3*a1+a1*a4+b3*b1+b1*b4+c3*c1+c1*c4)/rdet/360.0;
    C_G[6][0] = (a1*a1+a1*a4+b1*b1+b1*b4+c1*c1+c1*c4)/rdet/24.0;
    C_G[6][1] = (a2*a1+a2*a4+b2*b1+b2*b4+c2*c1+c2*c4)/rdet/24.0;
    C_G[6][2] = (a3*a1+a3*a4+b3*b1+b3*b4+c3*c1+c3*c4)/rdet/24.0;
    C_G[6][3] = (a4*a4+a1*a4+b4*b4+b1*b4+c4*c4+c1*c4)/rdet/24.0;
    C_G[6][4] = (a1*a4+a1*a1+2.0*a2*a4+a2*a1+b1*b4+b1*b1+2.0*b2*b4+b2*b1+c1*
        c4+c1*c1+2.0*c2*c4+c2*c1)/rdet/120.0;
    C_G[6][5] = (2.0*a3*a4+a3*a1+a1*a4+a1*a1+2.0*b3*b4+b3*b1+b1*b4+b1*b1+2.0*
        c3*c4+c3*c1+c1*c4+c1*c1)/rdet/120.0;
    C_G[6][6] = (a4*a4+a1*a4+a1*a1+b4*b4+b1*b4+b1*b1+c4*c4+c1*c4+c1*c1)/rdet/
        60.0;
    C_G[6][7] = (a2*a4+a3*a4+a2*a1+a3*a1+b2*b4+b3*b4+b2*b1+b3*b1+c2*c4+c3*c4+
        c2*c1+c3*c1)/rdet/120.0;
    C_G[6][8] = (a2*a4+a4*a4+2.0*a2*a1+a1*a4+b2*b4+b4*b4+2.0*b2*b1+b1*b4+c2*
        c4+c4*c4+2.0*c2*c1+c1*c4)/rdet/120.0;
    C_G[6][9] = (a3*a4+a4*a4+2.0*a3*a1+a1*a4+b3*b4+b4*b4+2.0*b3*b1+b1*b4+c3*
        c4+c4*c4+2.0*c3*c1+c1*c4)/rdet/120.0;
    C_G[6][10] = (2.0*a2*a4+2.0*a3*a4+a1*a4+a2*a1+a3*a1+a1*a1+2.0*b2*b4+2.0*
        b3*b4+b1*b4+b2*b1+b3*b1+b1*b1+2.0*c2*c4+2.0*c3*c4+c1*c4+c2*c1+c3*c1+c1*c1)/rdet
        /720.0;
    C_G[6][11] = (a2*a4+a4*a4+a1*a4+a2*a1+a1*a1+b2*b4+b4*b4+b1*b4+b2*b1+b1*b1
        +c2*c4+c4*c4+c1*c4+c2*c1+c1*c1)/rdet/360.0;
    C_G[6][12] = (a3*a4+a4*a4+a1*a4+a3*a1+a1*a1+b3*b4+b4*b4+b1*b4+b3*b1+b1*b1
        +c3*c4+c4*c4+c1*c4+c3*c1+c1*c1)/rdet/360.0;
    C_G[6][13] = (a3*a4+a4*a4+a2*a4+2.0*a3*a1+a1*a4+2.0*a2*a1+b3*b4+b4*b4+b2*
        b4+2.0*b3*b1+b1*b4+2.0*b2*b1+c3*c4+c4*c4+c2*c4+2.0*c3*c1+c1*c4+2.0*c2*c1)/rdet/
        720.0;
    C_G[6][14] = (a1*a4+a2*a4+b1*b4+b2*b4+c1*c4+c2*c4)/rdet/360.0;
    C_G[6][15] = (a3*a4+a1*a4+b3*b4+b1*b4+c3*c4+c1*c4)/rdet/360.0;
    C_G[6][16] = -(-a4*a4+a1*a1-b4*b4+b1*b1-c4*c4+c1*c1)/rdet/360.0;
    C_G[6][17] = 0.0;
    C_G[6][18] = -(a2*a1+a1*a4+b2*b1+b1*b4+c2*c1+c1*c4)/rdet/360.0;
    C_G[6][19] = -(a3*a1+a1*a4+b3*b1+b1*b4+c3*c1+c1*c4)/rdet/360.0;
    C_G[7][0] = (a2*a1+a3*a1+b2*b1+b3*b1+c2*c1+c3*c1)/rdet/24.0;
    C_G[7][1] = (a2*a2+a3*a2+b2*b2+b3*b2+c2*c2+c3*c2)/rdet/24.0;
    C_G[7][2] = (a3*a2+a3*a3+b3*b2+b3*b3+c3*c2+c3*c3)/rdet/24.0;
    C_G[7][3] = (a2*a4+a3*a4+b2*b4+b3*b4+c2*c4+c3*c4)/rdet/24.0;
    C_G[7][4] = (a2*a1+2.0*a3*a1+a2*a2+a3*a2+b2*b1+2.0*b3*b1+b2*b2+b3*b2+c2*
        c1+2.0*c3*c1+c2*c2+c3*c2)/rdet/120.0;
    C_G[7][5] = (a3*a2+a3*a3+2.0*a2*a1+a3*a1+b3*b2+b3*b3+2.0*b2*b1+b3*b1+c3*
        c2+c3*c3+2.0*c2*c1+c3*c1)/rdet/120.0;
    C_G[7][6] = (a2*a4+a3*a4+a2*a1+a3*a1+b2*b4+b3*b4+b2*b1+b3*b1+c2*c4+c3*c4+
        c2*c1+c3*c1)/rdet/120.0;
    C_G[7][7] = (a2*a2+a3*a2+a3*a3+b2*b2+b3*b2+b3*b3+c2*c2+c3*c2+c3*c3)/rdet/
        60.0;
    C_G[7][8] = (a2*a2+a2*a4+a3*a2+2.0*a3*a4+b2*b2+b2*b4+b3*b2+2.0*b3*b4+c2*
        c2+c2*c4+c3*c2+2.0*c3*c4)/rdet/120.0;
    C_G[7][9] = (a3*a2+2.0*a2*a4+a3*a3+a3*a4+b3*b2+2.0*b2*b4+b3*b3+b3*b4+c3*
        c2+2.0*c2*c4+c3*c3+c3*c4)/rdet/120.0;
    C_G[7][10] = (a2*a2+a3*a2+a2*a1+a3*a3+a3*a1+b2*b2+b3*b2+b2*b1+b3*b3+b3*b1
        +c2*c2+c3*c2+c2*c1+c3*c3+c3*c1)/rdet/360.0;
    C_G[7][11] = (a2*a2+a2*a4+a2*a1+a3*a2+2.0*a3*a4+2.0*a3*a1+b2*b2+b2*b4+b2*
        b1+b3*b2+2.0*b3*b4+2.0*b3*b1+c2*c2+c2*c4+c2*c1+c3*c2+2.0*c3*c4+2.0*c3*c1)/rdet/
        720.0;
    C_G[7][12] = (a3*a2+2.0*a2*a4+2.0*a2*a1+a3*a3+a3*a4+a3*a1+b3*b2+2.0*b2*b4
        +2.0*b2*b1+b3*b3+b3*b4+b3*b1+c3*c2+2.0*c2*c4+2.0*c2*c1+c3*c3+c3*c4+c3*c1)/rdet/
        720.0;
    C_G[7][13] = (a3*a2+a2*a4+a2*a2+a3*a3+a3*a4+b3*b2+b2*b4+b2*b2+b3*b3+b3*b4
        +c3*c2+c2*c4+c2*c2+c3*c3+c3*c4)/rdet/360.0;
    C_G[7][14] = -(a3*a1+a3*a2+b3*b1+b3*b2+c3*c1+c3*c2)/rdet/360.0;
    C_G[7][15] = -(a2*a1+a3*a2+b2*b1+b3*b2+c2*c1+c3*c2)/rdet/360.0;
    C_G[7][16] = 0.0;
    C_G[7][17] = -(a2*a2-a3*a3+b2*b2-b3*b3+c2*c2-c3*c3)/rdet/360.0;
    C_G[7][18] = (a3*a2+a3*a4+b3*b2+b3*b4+c3*c2+c3*c4)/rdet/360.0;
    C_G[7][19] = (a3*a2+a2*a4+b3*b2+b2*b4+c3*c2+c2*c4)/rdet/360.0;
    C_G[8][0] = (a2*a1+a1*a4+b2*b1+b1*b4+c2*c1+c1*c4)/rdet/24.0;
    C_G[8][1] = (a2*a2+a2*a4+b2*b2+b2*b4+c2*c2+c2*c4)/rdet/24.0;
    C_G[8][2] = (a3*a2+a3*a4+b3*b2+b3*b4+c3*c2+c3*c4)/rdet/24.0;
    C_G[8][3] = (a2*a4+a4*a4+b2*b4+b4*b4+c2*c4+c4*c4)/rdet/24.0;
    C_G[8][4] = (a2*a1+2.0*a1*a4+a2*a2+a2*a4+b2*b1+2.0*b1*b4+b2*b2+b2*b4+c2*
        c1+2.0*c1*c4+c2*c2+c2*c4)/rdet/120.0;
    C_G[8][5] = (a3*a2+a3*a4+a2*a1+a1*a4+b3*b2+b3*b4+b2*b1+b1*b4+c3*c2+c3*c4+
        c2*c1+c1*c4)/rdet/120.0;
    C_G[8][6] = (a2*a4+a4*a4+2.0*a2*a1+a1*a4+b2*b4+b4*b4+2.0*b2*b1+b1*b4+c2*
        c4+c4*c4+2.0*c2*c1+c1*c4)/rdet/120.0;
    C_G[8][7] = (a2*a2+a2*a4+a3*a2+2.0*a3*a4+b2*b2+b2*b4+b3*b2+2.0*b3*b4+c2*
        c2+c2*c4+c3*c2+2.0*c3*c4)/rdet/120.0;
    C_G[8][8] = (a2*a2+a2*a4+a4*a4+b2*b2+b2*b4+b4*b4+c2*c2+c2*c4+c4*c4)/rdet/
        60.0;
    C_G[8][9] = (2.0*a3*a2+a2*a4+a3*a4+a4*a4+2.0*b3*b2+b2*b4+b3*b4+b4*b4+2.0*
        c3*c2+c2*c4+c3*c4+c4*c4)/rdet/120.0;
    C_G[8][10] = (a2*a2+a3*a2+a2*a1+a2*a4+2.0*a3*a4+2.0*a1*a4+b2*b2+b3*b2+b2*
        b1+b2*b4+2.0*b3*b4+2.0*b1*b4+c2*c2+c3*c2+c2*c1+c2*c4+2.0*c3*c4+2.0*c1*c4)/rdet/
        720.0;
    C_G[8][11] = (a2*a2+a2*a4+a2*a1+a4*a4+a1*a4+b2*b2+b2*b4+b2*b1+b4*b4+b1*b4
        +c2*c2+c2*c4+c2*c1+c4*c4+c1*c4)/rdet/360.0;
    C_G[8][12] = (2.0*a3*a2+a2*a4+2.0*a2*a1+a3*a4+a4*a4+a1*a4+2.0*b3*b2+b2*b4
        +2.0*b2*b1+b3*b4+b4*b4+b1*b4+2.0*c3*c2+c2*c4+2.0*c2*c1+c3*c4+c4*c4+c1*c4)/rdet/
        720.0;
    C_G[8][13] = (a3*a2+a2*a4+a2*a2+a3*a4+a4*a4+b3*b2+b2*b4+b2*b2+b3*b4+b4*b4
        +c3*c2+c2*c4+c2*c2+c3*c4+c4*c4)/rdet/360.0;
    C_G[8][14] = -(a1*a4+a2*a4+b1*b4+b2*b4+c1*c4+c2*c4)/rdet/360.0;
    C_G[8][15] = 0.0;
    C_G[8][16] = -(a2*a1+a2*a4+b2*b1+b2*b4+c2*c1+c2*c4)/rdet/360.0;
    C_G[8][17] = (a2*a4+a3*a4+b2*b4+b3*b4+c2*c4+c3*c4)/rdet/360.0;
    C_G[8][18] = -(a2*a2-a4*a4+b2*b2-b4*b4+c2*c2-c4*c4)/rdet/360.0;
    C_G[8][19] = -(a3*a2+a2*a4+b3*b2+b2*b4+c3*c2+c2*c4)/rdet/360.0;
    C_G[9][0] = (a3*a1+a1*a4+b3*b1+b1*b4+c3*c1+c1*c4)/rdet/24.0;
    C_G[9][1] = (a3*a2+a2*a4+b3*b2+b2*b4+c3*c2+c2*c4)/rdet/24.0;
    C_G[9][2] = (a3*a3+a3*a4+b3*b3+b3*b4+c3*c3+c3*c4)/rdet/24.0;
    C_G[9][3] = (a3*a4+a4*a4+b3*b4+b4*b4+c3*c4+c4*c4)/rdet/24.0;
    C_G[9][4] = (a3*a1+a1*a4+a3*a2+a2*a4+b3*b1+b1*b4+b3*b2+b2*b4+c3*c1+c1*c4+
        c3*c2+c2*c4)/rdet/120.0;
    C_G[9][5] = (a3*a3+a3*a4+a3*a1+2.0*a1*a4+b3*b3+b3*b4+b3*b1+2.0*b1*b4+c3*
        c3+c3*c4+c3*c1+2.0*c1*c4)/rdet/120.0;
    C_G[9][6] = (a3*a4+a4*a4+2.0*a3*a1+a1*a4+b3*b4+b4*b4+2.0*b3*b1+b1*b4+c3*
        c4+c4*c4+2.0*c3*c1+c1*c4)/rdet/120.0;
    C_G[9][7] = (a3*a2+2.0*a2*a4+a3*a3+a3*a4+b3*b2+2.0*b2*b4+b3*b3+b3*b4+c3*
        c2+2.0*c2*c4+c3*c3+c3*c4)/rdet/120.0;
    C_G[9][8] = (2.0*a3*a2+a2*a4+a3*a4+a4*a4+2.0*b3*b2+b2*b4+b3*b4+b4*b4+2.0*
        c3*c2+c2*c4+c3*c4+c4*c4)/rdet/120.0;
    C_G[9][9] = (a3*a3+a3*a4+a4*a4+b3*b3+b3*b4+b4*b4+c3*c3+c3*c4+c4*c4)/rdet/
        60.0;
    C_G[9][10] = (a3*a2+a3*a3+a3*a1+2.0*a2*a4+a3*a4+2.0*a1*a4+b3*b2+b3*b3+b3*
        b1+2.0*b2*b4+b3*b4+2.0*b1*b4+c3*c2+c3*c3+c3*c1+2.0*c2*c4+c3*c4+2.0*c1*c4)/rdet/
        720.0;
    C_G[9][11] = (2.0*a3*a2+a3*a4+2.0*a3*a1+a2*a4+a4*a4+a1*a4+2.0*b3*b2+b3*b4
        +2.0*b3*b1+b2*b4+b4*b4+b1*b4+2.0*c3*c2+c3*c4+2.0*c3*c1+c2*c4+c4*c4+c1*c4)/rdet/
        720.0;
    C_G[9][12] = (a3*a3+a3*a4+a3*a1+a4*a4+a1*a4+b3*b3+b3*b4+b3*b1+b4*b4+b1*b4
        +c3*c3+c3*c4+c3*c1+c4*c4+c1*c4)/rdet/360.0;
    C_G[9][13] = (a3*a3+a3*a4+a3*a2+a4*a4+a2*a4+b3*b3+b3*b4+b3*b2+b4*b4+b2*b4
        +c3*c3+c3*c4+c3*c2+c4*c4+c2*c4)/rdet/360.0;
    C_G[9][14] = 0.0;
    C_G[9][15] = -(a3*a4+a1*a4+b3*b4+b1*b4+c3*c4+c1*c4)/rdet/360.0;
    C_G[9][16] = -(a3*a1+a3*a4+b3*b1+b3*b4+c3*c1+c3*c4)/rdet/360.0;
    C_G[9][17] = -(a2*a4+a3*a4+b2*b4+b3*b4+c2*c4+c3*c4)/rdet/360.0;
    C_G[9][18] = -(a3*a2+a3*a4+b3*b2+b3*b4+c3*c2+c3*c4)/rdet/360.0;
    C_G[9][19] = -(a3*a3-a4*a4+b3*b3-b4*b4+c3*c3-c4*c4)/rdet/360.0;
    C_G[10][0] = (a2*a1+a3*a1+a1*a1+b2*b1+b3*b1+b1*b1+c2*c1+c3*c1+c1*c1)/rdet
        /120.0;
    C_G[10][1] = (a2*a2+a3*a2+a2*a1+b2*b2+b3*b2+b2*b1+c2*c2+c3*c2+c2*c1)/rdet
        /120.0;
    C_G[10][2] = (a3*a2+a3*a3+a3*a1+b3*b2+b3*b3+b3*b1+c3*c2+c3*c3+c3*c1)/rdet
        /120.0;
    C_G[10][3] = (a2*a4+a3*a4+a1*a4+b2*b4+b3*b4+b1*b4+c2*c4+c3*c4+c1*c4)/rdet
        /120.0;
    C_G[10][4] = (a2*a1+a3*a1+a1*a1+a2*a2+a3*a2+b2*b1+b3*b1+b1*b1+b2*b2+b3*b2
        +c2*c1+c3*c1+c1*c1+c2*c2+c3*c2)/rdet/360.0;
    C_G[10][5] = (a3*a2+a3*a3+a3*a1+a2*a1+a1*a1+b3*b2+b3*b3+b3*b1+b2*b1+b1*b1
        +c3*c2+c3*c3+c3*c1+c2*c1+c1*c1)/rdet/360.0;
    C_G[10][6] = (2.0*a2*a4+2.0*a3*a4+a1*a4+a2*a1+a3*a1+a1*a1+2.0*b2*b4+2.0*
        b3*b4+b1*b4+b2*b1+b3*b1+b1*b1+2.0*c2*c4+2.0*c3*c4+c1*c4+c2*c1+c3*c1+c1*c1)/rdet
        /720.0;
    C_G[10][7] = (a2*a2+a3*a2+a2*a1+a3*a3+a3*a1+b2*b2+b3*b2+b2*b1+b3*b3+b3*b1
        +c2*c2+c3*c2+c2*c1+c3*c3+c3*c1)/rdet/360.0;
    C_G[10][8] = (a2*a2+a3*a2+a2*a1+a2*a4+2.0*a3*a4+2.0*a1*a4+b2*b2+b3*b2+b2*
        b1+b2*b4+2.0*b3*b4+2.0*b1*b4+c2*c2+c3*c2+c2*c1+c2*c4+2.0*c3*c4+2.0*c1*c4)/rdet/
        720.0;
    C_G[10][9] = (a3*a2+a3*a3+a3*a1+2.0*a2*a4+a3*a4+2.0*a1*a4+b3*b2+b3*b3+b3*
        b1+2.0*b2*b4+b3*b4+2.0*b1*b4+c3*c2+c3*c3+c3*c1+2.0*c2*c4+c3*c4+2.0*c1*c4)/rdet/
        720.0;
    C_G[10][10] = (a2*a2+a3*a2+a2*a1+a3*a3+a3*a1+a1*a1+b2*b2+b3*b2+b2*b1+b3*
        b3+b3*b1+b1*b1+c2*c2+c3*c2+c2*c1+c3*c3+c3*c1+c1*c1)/rdet/1260.0;
    C_G[10][11] = (a2*a2+a2*a4+a2*a1+a3*a2+2.0*a3*a4+a3*a1+a1*a4+a1*a1+b2*b2+
        b2*b4+b2*b1+b3*b2+2.0*b3*b4+b3*b1+b1*b4+b1*b1+c2*c2+c2*c4+c2*c1+c3*c2+2.0*c3*c4
        +c3*c1+c1*c4+c1*c1)/rdet/2520.0;
    C_G[10][12] = (a3*a2+2.0*a2*a4+a2*a1+a3*a3+a3*a4+a3*a1+a1*a4+a1*a1+b3*b2+
        2.0*b2*b4+b2*b1+b3*b3+b3*b4+b3*b1+b1*b4+b1*b1+c3*c2+2.0*c2*c4+c2*c1+c3*c3+c3*c4
        +c3*c1+c1*c4+c1*c1)/rdet/2520.0;
    C_G[10][13] = (a3*a2+a2*a4+a2*a2+a3*a3+a3*a4+a3*a1+2.0*a1*a4+a2*a1+b3*b2+
        b2*b4+b2*b2+b3*b3+b3*b4+b3*b1+2.0*b1*b4+b2*b1+c3*c2+c2*c4+c2*c2+c3*c3+c3*c4+c3*
        c1+2.0*c1*c4+c2*c1)/rdet/2520.0;
    C_G[10][14] = -(a3*a2+b3*b2+a1*a1+c3*c2-a2*a2-b2*b2-c2*c2-a3*a1+b1*b1+c1*
        c1-b3*b1-c3*c1)/rdet/2520.0;
    C_G[10][15] = -(a3*a2+b3*b2+a1*a1+c3*c2-b2*b1-a2*a1-c2*c1+b1*b1+c1*c1-a3*
        a3-b3*b3-c3*c3)/rdet/2520.0;
    C_G[10][16] = (a2*a4+a3*a4+a2*a1+a3*a1+b2*b4+b3*b4+b2*b1+b3*b1+c2*c4+c3*
        c4+c2*c1+c3*c1)/rdet/2520.0;
    C_G[10][17] = (b2*b1+a2*a1+c2*c1-a2*a2-b2*b2-c2*c2-a3*a1-b3*b1+a3*a3+b3*
        b3+c3*c3-c3*c1)/rdet/2520.0;
    C_G[10][18] = (a3*a2+a3*a4+a2*a1+a1*a4+b3*b2+b3*b4+b2*b1+b1*b4+c3*c2+c3*
        c4+c2*c1+c1*c4)/rdet/2520.0;
    C_G[10][19] = (a3*a1+a1*a4+a3*a2+a2*a4+b3*b1+b1*b4+b3*b2+b2*b4+c3*c1+c1*
        c4+c3*c2+c2*c4)/rdet/2520.0;
    C_G[11][0] = (a2*a1+a1*a4+a1*a1+b2*b1+b1*b4+b1*b1+c2*c1+c1*c4+c1*c1)/rdet
        /120.0;
    C_G[11][1] = (a2*a2+a2*a4+a2*a1+b2*b2+b2*b4+b2*b1+c2*c2+c2*c4+c2*c1)/rdet
        /120.0;
    C_G[11][2] = (a3*a2+a3*a4+a3*a1+b3*b2+b3*b4+b3*b1+c3*c2+c3*c4+c3*c1)/rdet
        /120.0;
    C_G[11][3] = (a2*a4+a4*a4+a1*a4+b2*b4+b4*b4+b1*b4+c2*c4+c4*c4+c1*c4)/rdet
        /120.0;
    C_G[11][4] = (a2*a1+a1*a4+a1*a1+a2*a2+a2*a4+b2*b1+b1*b4+b1*b1+b2*b2+b2*b4
        +c2*c1+c1*c4+c1*c1+c2*c2+c2*c4)/rdet/360.0;
    C_G[11][5] = (2.0*a3*a2+2.0*a3*a4+a3*a1+a2*a1+a1*a4+a1*a1+2.0*b3*b2+2.0*
        b3*b4+b3*b1+b2*b1+b1*b4+b1*b1+2.0*c3*c2+2.0*c3*c4+c3*c1+c2*c1+c1*c4+c1*c1)/rdet
        /720.0;
    C_G[11][6] = (a2*a4+a4*a4+a1*a4+a2*a1+a1*a1+b2*b4+b4*b4+b1*b4+b2*b1+b1*b1
        +c2*c4+c4*c4+c1*c4+c2*c1+c1*c1)/rdet/360.0;
    C_G[11][7] = (a2*a2+a2*a4+a2*a1+a3*a2+2.0*a3*a4+2.0*a3*a1+b2*b2+b2*b4+b2*
        b1+b3*b2+2.0*b3*b4+2.0*b3*b1+c2*c2+c2*c4+c2*c1+c3*c2+2.0*c3*c4+2.0*c3*c1)/rdet/
        720.0;
    C_G[11][8] = (a2*a2+a2*a4+a2*a1+a4*a4+a1*a4+b2*b2+b2*b4+b2*b1+b4*b4+b1*b4
        +c2*c2+c2*c4+c2*c1+c4*c4+c1*c4)/rdet/360.0;
    C_G[11][9] = (2.0*a3*a2+a3*a4+2.0*a3*a1+a2*a4+a4*a4+a1*a4+2.0*b3*b2+b3*b4
        +2.0*b3*b1+b2*b4+b4*b4+b1*b4+2.0*c3*c2+c3*c4+2.0*c3*c1+c2*c4+c4*c4+c1*c4)/rdet/
        720.0;
    C_G[11][10] = (a2*a2+a2*a4+a2*a1+a3*a2+2.0*a3*a4+a3*a1+a1*a4+a1*a1+b2*b2+
        b2*b4+b2*b1+b3*b2+2.0*b3*b4+b3*b1+b1*b4+b1*b1+c2*c2+c2*c4+c2*c1+c3*c2+2.0*c3*c4
        +c3*c1+c1*c4+c1*c1)/rdet/2520.0;
    C_G[11][11] = (a2*a2+a2*a4+a2*a1+a4*a4+a1*a4+a1*a1+b2*b2+b2*b4+b2*b1+b4*
        b4+b1*b4+b1*b1+c2*c2+c2*c4+c2*c1+c4*c4+c1*c4+c1*c1)/rdet/1260.0;
    C_G[11][12] = (2.0*a3*a2+a2*a4+a2*a1+a3*a4+a4*a4+a1*a4+a3*a1+a1*a1+2.0*b3
        *b2+b2*b4+b2*b1+b3*b4+b4*b4+b1*b4+b3*b1+b1*b1+2.0*c3*c2+c2*c4+c2*c1+c3*c4+c4*c4
        +c1*c4+c3*c1+c1*c1)/rdet/2520.0;
    C_G[11][13] = (a3*a2+a2*a4+a2*a2+a3*a4+a4*a4+2.0*a3*a1+a1*a4+a2*a1+b3*b2+
        b2*b4+b2*b2+b3*b4+b4*b4+2.0*b3*b1+b1*b4+b2*b1+c3*c2+c2*c4+c2*c2+c3*c4+c4*c4+2.0
        *c3*c1+c1*c4+c2*c1)/rdet/2520.0;
    C_G[11][14] = -(-b1*b4-c1*c4+a1*a1-a2*a2-b2*b2-c2*c2-a1*a4+b1*b1+c1*c1+a2
        *a4+b2*b4+c2*c4)/rdet/2520.0;
    C_G[11][15] = (a3*a2+a3*a4+a2*a1+a1*a4+b3*b2+b3*b4+b2*b1+b1*b4+c3*c2+c3*
        c4+c2*c1+c1*c4)/rdet/2520.0;
    C_G[11][16] = -(a1*a1-b2*b1-a2*a1-c2*c1-a4*a4-b4*b4-c4*c4+b1*b1+c1*c1+a2*
        a4+b2*b4+c2*c4)/rdet/2520.0;
    C_G[11][17] = (a2*a4+a3*a4+a2*a1+a3*a1+b2*b4+b3*b4+b2*b1+b3*b1+c2*c4+c3*
        c4+c2*c1+c3*c1)/rdet/2520.0;
    C_G[11][18] = (-b1*b4-c1*c4+b2*b1+a2*a1+c2*c1+a4*a4+b4*b4-a2*a2-b2*b2-c2*
        c2+c4*c4-a1*a4)/rdet/2520.0;
    C_G[11][19] = -(a3*a1+a1*a4+a3*a2+a2*a4+b3*b1+b1*b4+b3*b2+b2*b4+c3*c1+c1*
        c4+c3*c2+c2*c4)/rdet/2520.0;
    C_G[12][0] = (a3*a1+a1*a4+a1*a1+b3*b1+b1*b4+b1*b1+c3*c1+c1*c4+c1*c1)/rdet
        /120.0;
    C_G[12][1] = (a3*a2+a2*a4+a2*a1+b3*b2+b2*b4+b2*b1+c3*c2+c2*c4+c2*c1)/rdet
        /120.0;
    C_G[12][2] = (a3*a3+a3*a4+a3*a1+b3*b3+b3*b4+b3*b1+c3*c3+c3*c4+c3*c1)/rdet
        /120.0;
    C_G[12][3] = (a3*a4+a4*a4+a1*a4+b3*b4+b4*b4+b1*b4+c3*c4+c4*c4+c1*c4)/rdet
        /120.0;
    C_G[12][4] = (a3*a1+a1*a4+a1*a1+2.0*a3*a2+2.0*a2*a4+a2*a1+b3*b1+b1*b4+b1*
        b1+2.0*b3*b2+2.0*b2*b4+b2*b1+c3*c1+c1*c4+c1*c1+2.0*c3*c2+2.0*c2*c4+c2*c1)/rdet/
        720.0;
    C_G[12][5] = (a3*a3+a3*a4+a3*a1+a1*a4+a1*a1+b3*b3+b3*b4+b3*b1+b1*b4+b1*b1
        +c3*c3+c3*c4+c3*c1+c1*c4+c1*c1)/rdet/360.0;
    C_G[12][6] = (a3*a4+a4*a4+a1*a4+a3*a1+a1*a1+b3*b4+b4*b4+b1*b4+b3*b1+b1*b1
        +c3*c4+c4*c4+c1*c4+c3*c1+c1*c1)/rdet/360.0;
    C_G[12][7] = (a3*a2+2.0*a2*a4+2.0*a2*a1+a3*a3+a3*a4+a3*a1+b3*b2+2.0*b2*b4
        +2.0*b2*b1+b3*b3+b3*b4+b3*b1+c3*c2+2.0*c2*c4+2.0*c2*c1+c3*c3+c3*c4+c3*c1)/rdet/
        720.0;
    C_G[12][8] = (2.0*a3*a2+a2*a4+2.0*a2*a1+a3*a4+a4*a4+a1*a4+2.0*b3*b2+b2*b4
        +2.0*b2*b1+b3*b4+b4*b4+b1*b4+2.0*c3*c2+c2*c4+2.0*c2*c1+c3*c4+c4*c4+c1*c4)/rdet/
        720.0;
    C_G[12][9] = (a3*a3+a3*a4+a3*a1+a4*a4+a1*a4+b3*b3+b3*b4+b3*b1+b4*b4+b1*b4
        +c3*c3+c3*c4+c3*c1+c4*c4+c1*c4)/rdet/360.0;
    C_G[12][10] = (a3*a2+2.0*a2*a4+a2*a1+a3*a3+a3*a4+a3*a1+a1*a4+a1*a1+b3*b2+
        2.0*b2*b4+b2*b1+b3*b3+b3*b4+b3*b1+b1*b4+b1*b1+c3*c2+2.0*c2*c4+c2*c1+c3*c3+c3*c4
        +c3*c1+c1*c4+c1*c1)/rdet/2520.0;
    C_G[12][11] = (2.0*a3*a2+a2*a4+a2*a1+a3*a4+a4*a4+a1*a4+a3*a1+a1*a1+2.0*b3
        *b2+b2*b4+b2*b1+b3*b4+b4*b4+b1*b4+b3*b1+b1*b1+2.0*c3*c2+c2*c4+c2*c1+c3*c4+c4*c4
        +c1*c4+c3*c1+c1*c1)/rdet/2520.0;
    C_G[12][12] = (a3*a3+a3*a4+a3*a1+a4*a4+a1*a4+a1*a1+b3*b3+b3*b4+b3*b1+b4*
        b4+b1*b4+b1*b1+c3*c3+c3*c4+c3*c1+c4*c4+c1*c4+c1*c1)/rdet/1260.0;
    C_G[12][13] = (a3*a3+a3*a4+a3*a2+a4*a4+a2*a4+a3*a1+a1*a4+2.0*a2*a1+b3*b3+
        b3*b4+b3*b2+b4*b4+b2*b4+b3*b1+b1*b4+2.0*b2*b1+c3*c3+c3*c4+c3*c2+c4*c4+c2*c4+c3*
        c1+c1*c4+2.0*c2*c1)/rdet/2520.0;
    C_G[12][14] = (a3*a1+a1*a4+a3*a2+a2*a4+b3*b1+b1*b4+b3*b2+b2*b4+c3*c1+c1*
        c4+c3*c2+c2*c4)/rdet/2520.0;
    C_G[12][15] = -(-b1*b4-c1*c4+a3*a4+a1*a1+b3*b4+c3*c4-a1*a4+b1*b1+c1*c1-a3
        *a3-b3*b3-c3*c3)/rdet/2520.0;
    C_G[12][16] = -(a3*a4+a1*a1+b3*b4-a4*a4-b4*b4+c3*c4-c4*c4-a3*a1+b1*b1+c1*
        c1-b3*b1-c3*c1)/rdet/2520.0;
    C_G[12][17] = -(a2*a4+a3*a4+a2*a1+a3*a1+b2*b4+b3*b4+b2*b1+b3*b1+c2*c4+c3*
        c4+c2*c1+c3*c1)/rdet/2520.0;
    C_G[12][18] = -(a3*a2+a3*a4+a2*a1+a1*a4+b3*b2+b3*b4+b2*b1+b1*b4+c3*c2+c3*
        c4+c2*c1+c1*c4)/rdet/2520.0;
    C_G[12][19] = (-b1*b4-c1*c4+a4*a4+b4*b4+c4*c4-a1*a4+a3*a1+b3*b1-a3*a3-b3*
        b3-c3*c3+c3*c1)/rdet/2520.0;
    C_G[13][0] = (a3*a1+a1*a4+a2*a1+b3*b1+b1*b4+b2*b1+c3*c1+c1*c4+c2*c1)/rdet
        /120.0;
    C_G[13][1] = (a3*a2+a2*a4+a2*a2+b3*b2+b2*b4+b2*b2+c3*c2+c2*c4+c2*c2)/rdet
        /120.0;
    C_G[13][2] = (a3*a3+a3*a4+a3*a2+b3*b3+b3*b4+b3*b2+c3*c3+c3*c4+c3*c2)/rdet
        /120.0;
    C_G[13][3] = (a3*a4+a4*a4+a2*a4+b3*b4+b4*b4+b2*b4+c3*c4+c4*c4+c2*c4)/rdet
        /120.0;
    C_G[13][4] = (2.0*a3*a1+2.0*a1*a4+a2*a1+a3*a2+a2*a4+a2*a2+2.0*b3*b1+2.0*
        b1*b4+b2*b1+b3*b2+b2*b4+b2*b2+2.0*c3*c1+2.0*c1*c4+c2*c1+c3*c2+c2*c4+c2*c2)/rdet
        /720.0;
    C_G[13][5] = (a3*a3+a3*a4+a3*a2+a3*a1+2.0*a1*a4+2.0*a2*a1+b3*b3+b3*b4+b3*
        b2+b3*b1+2.0*b1*b4+2.0*b2*b1+c3*c3+c3*c4+c3*c2+c3*c1+2.0*c1*c4+2.0*c2*c1)/rdet/
        720.0;
    C_G[13][6] = (a3*a4+a4*a4+a2*a4+2.0*a3*a1+a1*a4+2.0*a2*a1+b3*b4+b4*b4+b2*
        b4+2.0*b3*b1+b1*b4+2.0*b2*b1+c3*c4+c4*c4+c2*c4+2.0*c3*c1+c1*c4+2.0*c2*c1)/rdet/
        720.0;
    C_G[13][7] = (a3*a2+a2*a4+a2*a2+a3*a3+a3*a4+b3*b2+b2*b4+b2*b2+b3*b3+b3*b4
        +c3*c2+c2*c4+c2*c2+c3*c3+c3*c4)/rdet/360.0;
    C_G[13][8] = (a3*a2+a2*a4+a2*a2+a3*a4+a4*a4+b3*b2+b2*b4+b2*b2+b3*b4+b4*b4
        +c3*c2+c2*c4+c2*c2+c3*c4+c4*c4)/rdet/360.0;
    C_G[13][9] = (a3*a3+a3*a4+a3*a2+a4*a4+a2*a4+b3*b3+b3*b4+b3*b2+b4*b4+b2*b4
        +c3*c3+c3*c4+c3*c2+c4*c4+c2*c4)/rdet/360.0;
    C_G[13][10] = (a3*a2+a2*a4+a2*a2+a3*a3+a3*a4+a3*a1+2.0*a1*a4+a2*a1+b3*b2+
        b2*b4+b2*b2+b3*b3+b3*b4+b3*b1+2.0*b1*b4+b2*b1+c3*c2+c2*c4+c2*c2+c3*c3+c3*c4+c3*
        c1+2.0*c1*c4+c2*c1)/rdet/2520.0;
    C_G[13][11] = (a3*a2+a2*a4+a2*a2+a3*a4+a4*a4+2.0*a3*a1+a1*a4+a2*a1+b3*b2+
        b2*b4+b2*b2+b3*b4+b4*b4+2.0*b3*b1+b1*b4+b2*b1+c3*c2+c2*c4+c2*c2+c3*c4+c4*c4+2.0
        *c3*c1+c1*c4+c2*c1)/rdet/2520.0;
    C_G[13][12] = (a3*a3+a3*a4+a3*a2+a4*a4+a2*a4+a3*a1+a1*a4+2.0*a2*a1+b3*b3+
        b3*b4+b3*b2+b4*b4+b2*b4+b3*b1+b1*b4+2.0*b2*b1+c3*c3+c3*c4+c3*c2+c4*c4+c2*c4+c3*
        c1+c1*c4+2.0*c2*c1)/rdet/2520.0;
    C_G[13][13] = (a3*a3+a3*a4+a3*a2+a4*a4+a2*a4+a2*a2+b3*b3+b3*b4+b3*b2+b4*
        b4+b2*b4+b2*b2+c3*c3+c3*c4+c3*c2+c4*c4+c2*c4+c2*c2)/rdet/1260.0;
    C_G[13][14] = -(a3*a1+a1*a4+a3*a2+a2*a4+b3*b1+b1*b4+b3*b2+b2*b4+c3*c1+c1*
        c4+c3*c2+c2*c4)/rdet/2520.0;
    C_G[13][15] = -(a3*a2+a3*a4+a2*a1+a1*a4+b3*b2+b3*b4+b2*b1+b1*b4+c3*c2+c3*
        c4+c2*c1+c1*c4)/rdet/2520.0;
    C_G[13][16] = -(a2*a4+a3*a4+a2*a1+a3*a1+b2*b4+b3*b4+b2*b1+b3*b1+c2*c4+c3*
        c4+c2*c1+c3*c1)/rdet/2520.0;
    C_G[13][17] = -(b2*b2-c3*c3+a3*a4+c3*c4+b3*b4-a2*a4-b2*b4-c2*c4+c2*c2-a3*
        a3-b3*b3+a2*a2)/rdet/2520.0;
    C_G[13][18] = -(b2*b2-b3*b2-a3*a2+a3*a4+c3*c4-c3*c2+b3*b4-a4*a4-b4*b4-c4*
        c4+c2*c2+a2*a2)/rdet/2520.0;
    C_G[13][19] = (b3*b2+a3*a2-c3*c3+c3*c2-a2*a4+a4*a4+b4*b4-b2*b4+c4*c4-c2*
        c4-a3*a3-b3*b3)/rdet/2520.0;
    C_G[14][0] = 0.0;
    C_G[14][1] = 0.0;
    C_G[14][2] = 0.0;
    C_G[14][3] = 0.0;
    C_G[14][4] = -(a1*a1-a2*a2+b1*b1-b2*b2+c1*c1-c2*c2)/rdet/360.0;
    C_G[14][5] = (a3*a1+a3*a2+b3*b1+b3*b2+c3*c1+c3*c2)/rdet/360.0;
    C_G[14][6] = (a1*a4+a2*a4+b1*b4+b2*b4+c1*c4+c2*c4)/rdet/360.0;
    C_G[14][7] = -(a3*a1+a3*a2+b3*b1+b3*b2+c3*c1+c3*c2)/rdet/360.0;
    C_G[14][8] = -(a1*a4+a2*a4+b1*b4+b2*b4+c1*c4+c2*c4)/rdet/360.0;
    C_G[14][9] = 0.0;
    C_G[14][10] = -(a3*a2+b3*b2+a1*a1+c3*c2-a2*a2-b2*b2-c2*c2-a3*a1+b1*b1+c1*
        c1-b3*b1-c3*c1)/rdet/2520.0;
    C_G[14][11] = -(-b1*b4-c1*c4+a1*a1-a2*a2-b2*b2-c2*c2-a1*a4+b1*b1+c1*c1+a2
        *a4+b2*b4+c2*c4)/rdet/2520.0;
    C_G[14][12] = (a3*a1+a1*a4+a3*a2+a2*a4+b3*b1+b1*b4+b3*b2+b2*b4+c3*c1+c1*
        c4+c3*c2+c2*c4)/rdet/2520.0;
    C_G[14][13] = -(a3*a1+a1*a4+a3*a2+a2*a4+b3*b1+b1*b4+b3*b2+b2*b4+c3*c1+c1*
        c4+c3*c2+c2*c4)/rdet/2520.0;
    C_G[14][14] = (2.0*b2*b2+2.0*a1*a1+c2*c1+b2*b1+a2*a1+2.0*c2*c2+2.0*b1*b1+
        2.0*c1*c1+2.0*a2*a2)/rdet/630.0;
    C_G[14][15] = (a3*a1+a1*a1+2.0*a3*a2+a2*a1+b3*b1+b1*b1+2.0*b3*b2+b2*b1+c3
        *c1+c1*c1+2.0*c3*c2+c2*c1)/rdet/1260.0;
    C_G[14][16] = (a1*a4+a1*a1+2.0*a2*a4+a2*a1+b1*b4+b1*b1+2.0*b2*b4+b2*b1+c1
        *c4+c1*c1+2.0*c2*c4+c2*c1)/rdet/1260.0;
    C_G[14][17] = -(a2*a1+2.0*a3*a1+a2*a2+a3*a2+b2*b1+2.0*b3*b1+b2*b2+b3*b2+
        c2*c1+2.0*c3*c1+c2*c2+c3*c2)/rdet/1260.0;
    C_G[14][18] = -(a2*a1+2.0*a1*a4+a2*a2+a2*a4+b2*b1+2.0*b1*b4+b2*b2+b2*b4+
        c2*c1+2.0*c1*c4+c2*c2+c2*c4)/rdet/1260.0;
    C_G[14][19] = 0.0;
    C_G[15][0] = 0.0;
    C_G[15][1] = 0.0;
    C_G[15][2] = 0.0;
    C_G[15][3] = 0.0;
    C_G[15][4] = (a2*a1+a3*a2+b2*b1+b3*b2+c2*c1+c3*c2)/rdet/360.0;
    C_G[15][5] = -(-a3*a3+a1*a1-b3*b3+b1*b1-c3*c3+c1*c1)/rdet/360.0;
    C_G[15][6] = (a3*a4+a1*a4+b3*b4+b1*b4+c3*c4+c1*c4)/rdet/360.0;
    C_G[15][7] = -(a2*a1+a3*a2+b2*b1+b3*b2+c2*c1+c3*c2)/rdet/360.0;
    C_G[15][8] = 0.0;
    C_G[15][9] = -(a3*a4+a1*a4+b3*b4+b1*b4+c3*c4+c1*c4)/rdet/360.0;
    C_G[15][10] = -(a3*a2+b3*b2+a1*a1+c3*c2-b2*b1-a2*a1-c2*c1+b1*b1+c1*c1-a3*
        a3-b3*b3-c3*c3)/rdet/2520.0;
    C_G[15][11] = (a3*a2+a3*a4+a2*a1+a1*a4+b3*b2+b3*b4+b2*b1+b1*b4+c3*c2+c3*
        c4+c2*c1+c1*c4)/rdet/2520.0;
    C_G[15][12] = -(-b1*b4-c1*c4+a3*a4+a1*a1+b3*b4+c3*c4-a1*a4+b1*b1+c1*c1-a3
        *a3-b3*b3-c3*c3)/rdet/2520.0;
    C_G[15][13] = -(a3*a2+a3*a4+a2*a1+a1*a4+b3*b2+b3*b4+b2*b1+b1*b4+c3*c2+c3*
        c4+c2*c1+c1*c4)/rdet/2520.0;
    C_G[15][14] = (a3*a1+a1*a1+2.0*a3*a2+a2*a1+b3*b1+b1*b1+2.0*b3*b2+b2*b1+c3
        *c1+c1*c1+2.0*c3*c2+c2*c1)/rdet/1260.0;
    C_G[15][15] = (2.0*c3*c3+2.0*a1*a1+c3*c1+a3*a1+2.0*b1*b1+2.0*c1*c1+b3*b1+
        2.0*a3*a3+2.0*b3*b3)/rdet/630.0;
    C_G[15][16] = (2.0*a3*a4+a3*a1+a1*a4+a1*a1+2.0*b3*b4+b3*b1+b1*b4+b1*b1+
        2.0*c3*c4+c3*c1+c1*c4+c1*c1)/rdet/1260.0;
    C_G[15][17] = (a3*a2+a3*a3+2.0*a2*a1+a3*a1+b3*b2+b3*b3+2.0*b2*b1+b3*b1+c3
        *c2+c3*c3+2.0*c2*c1+c3*c1)/rdet/1260.0;
    C_G[15][18] = 0.0;
    C_G[15][19] = -(a3*a3+a3*a4+a3*a1+2.0*a1*a4+b3*b3+b3*b4+b3*b1+2.0*b1*b4+
        c3*c3+c3*c4+c3*c1+2.0*c1*c4)/rdet/1260.0;
    C_G[16][0] = 0.0;
    C_G[16][1] = 0.0;
    C_G[16][2] = 0.0;
    C_G[16][3] = 0.0;
    C_G[16][4] = (a2*a1+a2*a4+b2*b1+b2*b4+c2*c1+c2*c4)/rdet/360.0;
    C_G[16][5] = (a3*a1+a3*a4+b3*b1+b3*b4+c3*c1+c3*c4)/rdet/360.0;
    C_G[16][6] = -(-a4*a4+a1*a1-b4*b4+b1*b1-c4*c4+c1*c1)/rdet/360.0;
    C_G[16][7] = 0.0;
    C_G[16][8] = -(a2*a1+a2*a4+b2*b1+b2*b4+c2*c1+c2*c4)/rdet/360.0;
    C_G[16][9] = -(a3*a1+a3*a4+b3*b1+b3*b4+c3*c1+c3*c4)/rdet/360.0;
    C_G[16][10] = (a2*a4+a3*a4+a2*a1+a3*a1+b2*b4+b3*b4+b2*b1+b3*b1+c2*c4+c3*
        c4+c2*c1+c3*c1)/rdet/2520.0;
    C_G[16][11] = -(a1*a1-b2*b1-a2*a1-c2*c1-a4*a4-b4*b4-c4*c4+b1*b1+c1*c1+a2*
        a4+b2*b4+c2*c4)/rdet/2520.0;
    C_G[16][12] = -(a3*a4+a1*a1+b3*b4-a4*a4-b4*b4+c3*c4-c4*c4-a3*a1+b1*b1+c1*
        c1-b3*b1-c3*c1)/rdet/2520.0;
    C_G[16][13] = -(a2*a4+a3*a4+a2*a1+a3*a1+b2*b4+b3*b4+b2*b1+b3*b1+c2*c4+c3*
        c4+c2*c1+c3*c1)/rdet/2520.0;
    C_G[16][14] = (a1*a4+a1*a1+2.0*a2*a4+a2*a1+b1*b4+b1*b1+2.0*b2*b4+b2*b1+c1
        *c4+c1*c1+2.0*c2*c4+c2*c1)/rdet/1260.0;
    C_G[16][15] = (2.0*a3*a4+a3*a1+a1*a4+a1*a1+2.0*b3*b4+b3*b1+b1*b4+b1*b1+
        2.0*c3*c4+c3*c1+c1*c4+c1*c1)/rdet/1260.0;
    C_G[16][16] = (b1*b4+2.0*a1*a1+c1*c4+2.0*a4*a4+2.0*b4*b4+a1*a4+2.0*c4*c4+
        2.0*b1*b1+2.0*c1*c1)/rdet/630.0;
    C_G[16][17] = 0.0;
    C_G[16][18] = (a2*a4+a4*a4+2.0*a2*a1+a1*a4+b2*b4+b4*b4+2.0*b2*b1+b1*b4+c2
        *c4+c4*c4+2.0*c2*c1+c1*c4)/rdet/1260.0;
    C_G[16][19] = (a3*a4+a4*a4+2.0*a3*a1+a1*a4+b3*b4+b4*b4+2.0*b3*b1+b1*b4+c3
        *c4+c4*c4+2.0*c3*c1+c1*c4)/rdet/1260.0;
    C_G[17][0] = 0.0;
    C_G[17][1] = 0.0;
    C_G[17][2] = 0.0;
    C_G[17][3] = 0.0;
    C_G[17][4] = (a2*a1+a3*a1+b2*b1+b3*b1+c2*c1+c3*c1)/rdet/360.0;
    C_G[17][5] = -(a2*a1+a3*a1+b2*b1+b3*b1+c2*c1+c3*c1)/rdet/360.0;
    C_G[17][6] = 0.0;
    C_G[17][7] = -(a2*a2-a3*a3+b2*b2-b3*b3+c2*c2-c3*c3)/rdet/360.0;
    C_G[17][8] = (a2*a4+a3*a4+b2*b4+b3*b4+c2*c4+c3*c4)/rdet/360.0;
    C_G[17][9] = -(a2*a4+a3*a4+b2*b4+b3*b4+c2*c4+c3*c4)/rdet/360.0;
    C_G[17][10] = (b2*b1+a2*a1+c2*c1-a2*a2-b2*b2-c2*c2-a3*a1-b3*b1+a3*a3+b3*
        b3+c3*c3-c3*c1)/rdet/2520.0;
    C_G[17][11] = (a2*a4+a3*a4+a2*a1+a3*a1+b2*b4+b3*b4+b2*b1+b3*b1+c2*c4+c3*
        c4+c2*c1+c3*c1)/rdet/2520.0;
    C_G[17][12] = -(a2*a4+a3*a4+a2*a1+a3*a1+b2*b4+b3*b4+b2*b1+b3*b1+c2*c4+c3*
        c4+c2*c1+c3*c1)/rdet/2520.0;
    C_G[17][13] = -(b2*b2-c3*c3+a3*a4+c3*c4+b3*b4-a2*a4-b2*b4-c2*c4+c2*c2-a3*
        a3-b3*b3+a2*a2)/rdet/2520.0;
    C_G[17][14] = -(a2*a1+2.0*a3*a1+a2*a2+a3*a2+b2*b1+2.0*b3*b1+b2*b2+b3*b2+
        c2*c1+2.0*c3*c1+c2*c2+c3*c2)/rdet/1260.0;
    C_G[17][15] = (a3*a2+a3*a3+2.0*a2*a1+a3*a1+b3*b2+b3*b3+2.0*b2*b1+b3*b1+c3
        *c2+c3*c3+2.0*c2*c1+c3*c1)/rdet/1260.0;
    C_G[17][16] = 0.0;
    C_G[17][17] = (2.0*b3*b3+a3*a2+2.0*a2*a2+2.0*c3*c3+b3*b2+2.0*b2*b2+c3*c2+
        2.0*a3*a3+2.0*c2*c2)/rdet/630.0;
    C_G[17][18] = (a2*a2+a2*a4+a3*a2+2.0*a3*a4+b2*b2+b2*b4+b3*b2+2.0*b3*b4+c2
        *c2+c2*c4+c3*c2+2.0*c3*c4)/rdet/1260.0;
    C_G[17][19] = -(a3*a2+2.0*a2*a4+a3*a3+a3*a4+b3*b2+2.0*b2*b4+b3*b3+b3*b4+
        c3*c2+2.0*c2*c4+c3*c3+c3*c4)/rdet/1260.0;
    C_G[18][0] = 0.0;
    C_G[18][1] = 0.0;
    C_G[18][2] = 0.0;
    C_G[18][3] = 0.0;
    C_G[18][4] = (a2*a1+a1*a4+b2*b1+b1*b4+c2*c1+c1*c4)/rdet/360.0;
    C_G[18][5] = 0.0;
    C_G[18][6] = -(a2*a1+a1*a4+b2*b1+b1*b4+c2*c1+c1*c4)/rdet/360.0;
    C_G[18][7] = (a3*a2+a3*a4+b3*b2+b3*b4+c3*c2+c3*c4)/rdet/360.0;
    C_G[18][8] = -(a2*a2-a4*a4+b2*b2-b4*b4+c2*c2-c4*c4)/rdet/360.0;
    C_G[18][9] = -(a3*a2+a3*a4+b3*b2+b3*b4+c3*c2+c3*c4)/rdet/360.0;
    C_G[18][10] = (a3*a2+a3*a4+a2*a1+a1*a4+b3*b2+b3*b4+b2*b1+b1*b4+c3*c2+c3*
        c4+c2*c1+c1*c4)/rdet/2520.0;
    C_G[18][11] = (-b1*b4-c1*c4+b2*b1+a2*a1+c2*c1+a4*a4+b4*b4-a2*a2-b2*b2-c2*
        c2+c4*c4-a1*a4)/rdet/2520.0;
    C_G[18][12] = -(a3*a2+a3*a4+a2*a1+a1*a4+b3*b2+b3*b4+b2*b1+b1*b4+c3*c2+c3*
        c4+c2*c1+c1*c4)/rdet/2520.0;
    C_G[18][13] = -(b2*b2-b3*b2-a3*a2+a3*a4+c3*c4-c3*c2+b3*b4-a4*a4-b4*b4-c4*
        c4+c2*c2+a2*a2)/rdet/2520.0;
    C_G[18][14] = -(a2*a1+2.0*a1*a4+a2*a2+a2*a4+b2*b1+2.0*b1*b4+b2*b2+b2*b4+
        c2*c1+2.0*c1*c4+c2*c2+c2*c4)/rdet/1260.0;
    C_G[18][15] = 0.0;
    C_G[18][16] = (a2*a4+a4*a4+2.0*a2*a1+a1*a4+b2*b4+b4*b4+2.0*b2*b1+b1*b4+c2
        *c4+c4*c4+2.0*c2*c1+c1*c4)/rdet/1260.0;
    C_G[18][17] = (a2*a2+a2*a4+a3*a2+2.0*a3*a4+b2*b2+b2*b4+b3*b2+2.0*b3*b4+c2
        *c2+c2*c4+c3*c2+2.0*c3*c4)/rdet/1260.0;
    C_G[18][18] = (2.0*a4*a4+2.0*b4*b4+2.0*c4*c4+a2*a4+2.0*a2*a2+2.0*b2*b2+b2
        *b4+c2*c4+2.0*c2*c2)/rdet/630.0;
    C_G[18][19] = (2.0*a3*a2+a2*a4+a3*a4+a4*a4+2.0*b3*b2+b2*b4+b3*b4+b4*b4+
        2.0*c3*c2+c2*c4+c3*c4+c4*c4)/rdet/1260.0;
    C_G[19][0] = 0.0;
    C_G[19][1] = 0.0;
    C_G[19][2] = 0.0;
    C_G[19][3] = 0.0;
    C_G[19][4] = 0.0;
    C_G[19][5] = (a3*a1+a1*a4+b3*b1+b1*b4+c3*c1+c1*c4)/rdet/360.0;
    C_G[19][6] = -(a3*a1+a1*a4+b3*b1+b1*b4+c3*c1+c1*c4)/rdet/360.0;
    C_G[19][7] = (a3*a2+a2*a4+b3*b2+b2*b4+c3*c2+c2*c4)/rdet/360.0;
    C_G[19][8] = -(a3*a2+a2*a4+b3*b2+b2*b4+c3*c2+c2*c4)/rdet/360.0;
    C_G[19][9] = -(a3*a3-a4*a4+b3*b3-b4*b4+c3*c3-c4*c4)/rdet/360.0;
    C_G[19][10] = (a3*a1+a1*a4+a3*a2+a2*a4+b3*b1+b1*b4+b3*b2+b2*b4+c3*c1+c1*
        c4+c3*c2+c2*c4)/rdet/2520.0;
    C_G[19][11] = -(a3*a1+a1*a4+a3*a2+a2*a4+b3*b1+b1*b4+b3*b2+b2*b4+c3*c1+c1*
        c4+c3*c2+c2*c4)/rdet/2520.0;
    C_G[19][12] = (-b1*b4-c1*c4+a4*a4+b4*b4+c4*c4-a1*a4+a3*a1+b3*b1-a3*a3-b3*
        b3-c3*c3+c3*c1)/rdet/2520.0;
    C_G[19][13] = (b3*b2+a3*a2-c3*c3+c3*c2-a2*a4+a4*a4+b4*b4-b2*b4+c4*c4-c2*
        c4-a3*a3-b3*b3)/rdet/2520.0;
    C_G[19][14] = 0.0;
    C_G[19][15] = -(a3*a3+a3*a4+a3*a1+2.0*a1*a4+b3*b3+b3*b4+b3*b1+2.0*b1*b4+
        c3*c3+c3*c4+c3*c1+2.0*c1*c4)/rdet/1260.0;
    C_G[19][16] = (a3*a4+a4*a4+2.0*a3*a1+a1*a4+b3*b4+b4*b4+2.0*b3*b1+b1*b4+c3
        *c4+c4*c4+2.0*c3*c1+c1*c4)/rdet/1260.0;
    C_G[19][17] = -(a3*a2+2.0*a2*a4+a3*a3+a3*a4+b3*b2+2.0*b2*b4+b3*b3+b3*b4+
        c3*c2+2.0*c2*c4+c3*c3+c3*c4)/rdet/1260.0;
    C_G[19][18] = (2.0*a3*a2+a2*a4+a3*a4+a4*a4+2.0*b3*b2+b2*b4+b3*b4+b4*b4+
        2.0*c3*c2+c2*c4+c3*c4+c4*c4)/rdet/1260.0;
    C_G[19][19] = (a3*a4+b3*b4+2.0*a4*a4+2.0*b4*b4+2.0*c4*c4+2.0*b3*b3+c3*c4+
        2.0*c3*c3+2.0*a3*a3)/rdet/630.0;
/*
    int numRows = GetNodalBasisDim ();
    CRect_Matrix locStiff (numRows, numRows);

    for (int i = 0; i < numRows; i++) {
        for (int j = 0; j < numRows; j++) {
            locStiff (i, j) = C_G[i][j];
        }
    }
    return locStiff;
*/
    matrix_t<double, 10, 12> matr;
    for(size_t i = 0; i < 10; i++)
        for(size_t j = 0; j < 12; j++)
            matr[i][j] = C_G[i][j];
    return matr;
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
