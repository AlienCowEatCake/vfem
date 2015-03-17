#include "tetrahedron.h"

// ============================================================================

tetrahedron_base::tetrahedron_base()
{
    for(size_t i = 0; i < 4; i++)
        nodes[i] = NULL;
    for(size_t i = 0; i < 6; i++)
        edges[i] = NULL;
    phys = NULL;
}

const point & tetrahedron_base::get_node(size_t i) const
{
    if(!nodes[i])
    {
        cerr << "Error: Null pointer at get_node(" << i << ")" << endl;
        throw NULL_PTR_ERROR;
    }
    return (* nodes[i]);
}

const edge & tetrahedron_base::get_edge(size_t i) const
{
    if(!edges[i])
    {
        cerr << "Error: Null pointer at get_edge(" << i << ")" << endl;
        throw NULL_PTR_ERROR;
    }
    return (* edges[i]);
}

const phys_area & tetrahedron_base::get_phys_area() const
{
    if(!phys)
    {
        cerr << "Error: Null pointer at get_phys_area()" << endl;
        throw NULL_PTR_ERROR;
    }
    return * phys;
}

double tetrahedron_base::lambda(size_t i, const point & p) const
{
    return L[i][3] + L[i][0] * p.x + L[i][1] * p.y + L[i][2] * p.z;
}

vector3 tetrahedron_base::grad_lambda(size_t i) const
{
    return vector3(L[i][0], L[i][1], L[i][2]);
}

vector3 tetrahedron_base::w(size_t i, const point & p) const
{
    switch(i + 1)
    {
    case 1:
        return lambda(0, p) * grad_lambda(1) - lambda(1, p) * grad_lambda(0);
    case 2:
        return lambda(0, p) * grad_lambda(2) - lambda(2, p) * grad_lambda(0);
    case 3:
        return lambda(0, p) * grad_lambda(3) - lambda(3, p) * grad_lambda(0);
    case 4:
        return lambda(1, p) * grad_lambda(2) - lambda(2, p) * grad_lambda(1);
    case 5:
        return lambda(1, p) * grad_lambda(3) - lambda(3, p) * grad_lambda(1);
    case 6:
        return lambda(2, p) * grad_lambda(3) - lambda(3, p) * grad_lambda(2);
    case 7:
        return lambda(0, p) * grad_lambda(1) + lambda(1, p) * grad_lambda(0);
    case 8:
        return lambda(0, p) * grad_lambda(2) + lambda(2, p) * grad_lambda(0);
    case 9:
        return lambda(0, p) * grad_lambda(3) + lambda(3, p) * grad_lambda(0);
    case 10:
        return lambda(1, p) * grad_lambda(2) + lambda(2, p) * grad_lambda(1);
    case 11:
        return lambda(1, p) * grad_lambda(3) + lambda(3, p) * grad_lambda(1);
    case 12:
        return lambda(2, p) * grad_lambda(3) + lambda(3, p) * grad_lambda(2);
    }
    cerr << "Error: Incorrect basis function number!" << endl;
    throw(ADDRESSING_ERROR);
}

vector3 tetrahedron_base::rotw(size_t i, const point & p) const
{
    MAYBE_UNUSED(p);
    vector3 grad1, grad2;
    switch(i + 1)
    {
    case 1:
        grad1 = grad_lambda(0);
        grad2 = grad_lambda(1);
        break;
    case 2:
        grad1 = grad_lambda(0);
        grad2 = grad_lambda(2);
        break;
    case 3:
        grad1 = grad_lambda(0);
        grad2 = grad_lambda(3);
        break;
    case 4:
        grad1 = grad_lambda(1);
        grad2 = grad_lambda(2);
        break;
    case 5:
        grad1 = grad_lambda(1);
        grad2 = grad_lambda(3);
        break;
    case 6:
        grad1 = grad_lambda(2);
        grad2 = grad_lambda(3);
        break;
    case 7:
    case 8:
    case 9:
    case 10:
    case 11:
    case 12:
        return vector3(0.0, 0.0, 0.0);
    default:
        cerr << "Error: Incorrect rot basis function number!" << endl;
        throw(ADDRESSING_ERROR);
    }
    return 2.0 * grad1.cross(grad2);
}

void tetrahedron_base::init()
{
    matrix4 D;
    for(size_t i = 0; i < 3; i++)
        for(size_t j = 0; j < 4; j++)
            D[i][j] = get_node(j)[i];
    for(size_t j = 0; j < 4; j++)
        D[3][j] = 1.0;
    double D_det;
    L = inverse(D, D_det);

    jacobian = fabs(D_det);

    // Точки Гаусса на мастер-элементе
    static const double gauss_a = (5.0 - sqrt(5.0)) / 20.0;
    static const double gauss_b = (5.0 + 3.0 * sqrt(5.0)) / 20.0;
    double Gauss_cord[4][4];
    double Gauss_cord_gl[4][4];
    // Заполнение матрицы в локальных координатах
    Gauss_cord[0][0] = 1.0 - gauss_b - 2.0 * gauss_a;
    Gauss_cord[1][0] = gauss_b;
    Gauss_cord[2][0] = Gauss_cord[3][0] = gauss_a;
    Gauss_cord[0][1] = 1.0 - gauss_b - 2.0 * gauss_a;
    Gauss_cord[1][1] = Gauss_cord[3][1] = gauss_a;
    Gauss_cord[2][1] = gauss_b;
    Gauss_cord[0][2] = 1.0 - gauss_b - 2.0 * gauss_a;
    Gauss_cord[1][2] = Gauss_cord[2][2] = gauss_a;
    Gauss_cord[3][2] = gauss_b;
    Gauss_cord[0][3] = 1.0 - 3.0 * gauss_a;
    Gauss_cord[1][3] = Gauss_cord[2][3] = Gauss_cord[3][3] = gauss_a;

    // Перевод на текущий тетраэдр
    for(size_t i = 0; i < 4; i++)
    {
        for(size_t j = 0 ; j < 4; j++)
        {
            Gauss_cord_gl[i][j] = 0;
            for(size_t k = 0; k < 4; k++)
                Gauss_cord_gl[i][j] += D[i][k] * Gauss_cord[k][j];
        }
    }
    for(size_t i = 0; i < 4; i++)
        for(size_t j = 0; j < 3; j++)
            gauss_points[i][j] = Gauss_cord_gl[j][i];
    for(size_t i = 0; i < 4; i++)
        gauss_weights[i] = 1.0 / 24.0;

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
    double t[6][6];     //параметры прямых для всех прямых (6) и всех плоскостей (6)
    double cords_t[6][6][3];    //прямая, плоскость, координата
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

    for(size_t i = 0; i < 6; i++)   // берём прямоую и проверяем, что пересечение с плоскостью попадает в раcсматриваемый отрезок прямой и в рассатриваемую часть плоскости
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

double tetrahedron_base::diff_normL2(const carray12 & q, cvector3(*func)(const point &)) const
{
    complex<double> result = 0.0;
    for(size_t k = 0; k < 4; k++)
    {
        cvector3 val(0.0, 0.0, 0.0);
        for(size_t i = 0; i < 12; i++)
            val = val + q[i] * cvector3(w(i, gauss_points[k]));
        cvector3 func_d = func(gauss_points[k]) - val;
        cvector3 func_d_conj = func_d.cj();
        result += gauss_weights[k] * (func_d * func_d_conj);
    }
    result *= jacobian;
    return result.real();
}

double tetrahedron_base::diff_normL2(const carray12 & q, const carray12 & q_true) const
{
    complex<double> result = 0.0;
    for(size_t k = 0; k < 4; k++)
    {
        cvector3 val(0.0, 0.0, 0.0), val_true(0.0, 0.0, 0.0);
        for(size_t i = 0; i < 12; i++)
        {
            val = val + q[i] * cvector3(w(i, gauss_points[k]));
            val_true = val_true + q_true[i] * cvector3(w(i, gauss_points[k]));
        }
        cvector3 func_d = val_true - val;
        cvector3 func_d_conj = func_d.cj();
        result += gauss_weights[k] * (func_d * func_d_conj);
    }
    result *= jacobian;
    return result.real();
}

double tetrahedron_base::normL2(cvector3(*func)(const point &)) const
{
    complex<double> result = 0.0;
    for(size_t k = 0; k < 4; k++)
    {
        cvector3 func_d = func(gauss_points[k]);
        cvector3 func_d_conj = func_d.cj();
        result += gauss_weights[k] * (func_d * func_d_conj);
    }
    result *= jacobian;
    return result.real();
}

double tetrahedron_base::normL2(const carray12 & q_true) const
{
    complex<double> result = 0.0;
    for(size_t k = 0; k < 4; k++)
    {
        cvector3 func_d;
        for(size_t i = 0; i < 12; i++)
            func_d = func_d + q_true[i] * cvector3(w(i, gauss_points[k]));
        cvector3 func_d_conj = func_d.cj();
        result += gauss_weights[k] * (func_d * func_d_conj);
    }
    result *= jacobian;
    return result.real();
}

// ============================================================================

double tetrahedron::integrate_w(size_t i, size_t j) const
{
    double result = 0.0;
    for(size_t k = 0; k < 4; k++)
        result += gauss_weights[k] *
                  w(i, gauss_points[k]) *
                  w(j, gauss_points[k]);
    return result * jacobian;
}

double tetrahedron::integrate_rotw(size_t i, size_t j) const
{
    double result = 0.0;
    for(size_t k = 0; k < 4; k++)
        result += gauss_weights[k] *
                  rotw(i, gauss_points[k]) *
                  rotw(j, gauss_points[k]);
    return result * jacobian;
}

complex<double> tetrahedron::integrate_fw(cvector3(*func)(const point &, const phys_area &), size_t i) const
{
    complex<double> result = 0.0;
    for(size_t k = 0; k < 4; k++)
        result += gauss_weights[k] *
                  func(gauss_points[k], get_phys_area()) *
                  w(i, gauss_points[k]);
    return result * jacobian;
}

matrix12 tetrahedron::G() const
{
    matrix12 matr;
    for(size_t i = 0; i < 12; i++)
        for(size_t j = 0; j < 12; j++)
            matr[i][j] = 0.0;

    for(size_t i = 0; i < 6; i++)
        for(size_t j = 0; j <= i; j++)
            matr[j][i] = matr[i][j] = integrate_rotw(i, j);
    return matr;
}

matrix12 tetrahedron::M() const
{
    matrix12 matr;
    for(size_t i = 0; i < 12; i++)
        for(size_t j = 0; j <= i; j++)
            matr[j][i] = matr[i][j] = integrate_w(i, j);
    return matr;
}

carray12 tetrahedron::rp(cvector3(*func)(const point &, const phys_area &)) const
{
    carray12 arr;
    for(size_t i = 0; i < 12; i++)
        arr[i] = integrate_fw(func, i);
    return arr;
}

// ============================================================================
