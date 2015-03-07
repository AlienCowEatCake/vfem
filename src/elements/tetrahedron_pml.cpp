#include "tetrahedron.h"

const cpoint & tetrahedron_pml::get_node_pml(size_t i) const
{
    if(!nodes_pml[i])
    {
        cerr << "Error: Null pointer at get_node_pml(" << i << ")" << endl;
        throw NULL_PTR_ERROR;
    }
    return (* nodes_pml[i]);
}

void tetrahedron_pml::init_pml(cvector3(* get_s)(const point &, const tetrahedron_pml *, const phys_pml_area *), const phys_pml_area * phys_pml)
{
    this->get_s = get_s;
    this->phys_pml = phys_pml;

    cmatrix4 D;
    for(size_t i = 0; i < 3; i++)
        for(size_t j = 0; j < 4; j++)
            D[i][j] = get_node_pml(j)[i];
    for(size_t j = 0; j < 4; j++)
        D[3][j] = 1.0;
    complex<double> D_det;
    L_pml = inverse(D, D_det);

    jacobian_pml = abs(D_det);

    // Перевод точек Гаусса с мастер-элемента на текущий тетраэдр
    for(size_t i = 0; i < 3; i++)
    {
        for(size_t j = 0 ; j < tet_integration::gauss_num; j++)
        {
            gauss_points_pml[j][i] = 0;
            for(size_t k = 0; k < 4; k++)
                gauss_points_pml[j][i] += D[i][k] * tet_integration::gauss_points[j][k];
        }
    }
}

complex<double> tetrahedron_pml::lambda_pml(size_t i, const cpoint & p) const
{
    return L_pml[i][3] + L_pml[i][0] * p.x + L_pml[i][1] * p.y + L_pml[i][2] * p.z;
}

cvector3 tetrahedron_pml::grad_lambda_pml(size_t i) const
{
    return cvector3(L_pml[i][0], L_pml[i][1], L_pml[i][2]);
}

cvector3 tetrahedron_pml::w_pml(size_t i, const cpoint & p) const
{
    switch(i + 1)
    {
    case 1:
        return lambda_pml(0, p) * grad_lambda_pml(1) - lambda_pml(1, p) * grad_lambda_pml(0);
    case 2:
        return lambda_pml(0, p) * grad_lambda_pml(2) - lambda_pml(2, p) * grad_lambda_pml(0);
    case 3:
        return lambda_pml(0, p) * grad_lambda_pml(3) - lambda_pml(3, p) * grad_lambda_pml(0);
    case 4:
        return lambda_pml(1, p) * grad_lambda_pml(2) - lambda_pml(2, p) * grad_lambda_pml(1);
    case 5:
        return lambda_pml(1, p) * grad_lambda_pml(3) - lambda_pml(3, p) * grad_lambda_pml(1);
    case 6:
        return lambda_pml(2, p) * grad_lambda_pml(3) - lambda_pml(3, p) * grad_lambda_pml(2);
    case 7:
        return lambda_pml(0, p) * grad_lambda_pml(1) + lambda_pml(1, p) * grad_lambda_pml(0);
    case 8:
        return lambda_pml(0, p) * grad_lambda_pml(2) + lambda_pml(2, p) * grad_lambda_pml(0);
    case 9:
        return lambda_pml(0, p) * grad_lambda_pml(3) + lambda_pml(3, p) * grad_lambda_pml(0);
    case 10:
        return lambda_pml(1, p) * grad_lambda_pml(2) + lambda_pml(2, p) * grad_lambda_pml(1);
    case 11:
        return lambda_pml(1, p) * grad_lambda_pml(3) + lambda_pml(3, p) * grad_lambda_pml(1);
    case 12:
        return lambda_pml(2, p) * grad_lambda_pml(3) + lambda_pml(3, p) * grad_lambda_pml(2);
    }
    cerr << "Error: Incorrect basis function number!" << endl;
    throw(ADDRESSING_ERROR);
}

cvector3 tetrahedron_pml::rotw_pml(size_t i, const cpoint & p, const point & p_non_PML) const
{
    MAYBE_UNUSED(p);
    cvector3 grad1, grad2;
    switch(i + 1)
    {
    case 1:
        grad1 = grad_lambda_pml(0);
        grad2 = grad_lambda_pml(1);
        break;
    case 2:
        grad1 = grad_lambda_pml(0);
        grad2 = grad_lambda_pml(2);
        break;
    case 3:
        grad1 = grad_lambda_pml(0);
        grad2 = grad_lambda_pml(3);
        break;
    case 4:
        grad1 = grad_lambda_pml(1);
        grad2 = grad_lambda_pml(2);
        break;
    case 5:
        grad1 = grad_lambda_pml(1);
        grad2 = grad_lambda_pml(3);
        break;
    case 6:
        grad1 = grad_lambda_pml(2);
        grad2 = grad_lambda_pml(3);
        break;
    case 7:
    case 8:
    case 9:
    case 10:
    case 11:
    case 12:
        return cvector3(0.0, 0.0, 0.0);
    default:
        cerr << "Error: Incorrect rot basis function number!" << endl;
        throw(ADDRESSING_ERROR);
    }

    cvector3 s = get_s(p_non_PML, this, phys_pml);
    for(size_t k = 0; k < 3; k++)
    {
        grad1[k] /= s[k];
        grad2[k] /= s[k];
    }

    return 2.0 * grad1.cross(grad2);
}

complex<double> tetrahedron_pml::integrate_w(size_t i, size_t j) const
{
//#if defined __GNUC__
//#warning conjugate
//#endif
    complex<double> result = 0.0;
    for(size_t k = 0; k < tet_integration::gauss_num; k++)
        result += tet_integration::gauss_weights[k] *
                  w_pml(i, gauss_points_pml[k]).cj() *
                  w_pml(j, gauss_points_pml[k]).cj();
    return result * jacobian_pml;
}

complex<double> tetrahedron_pml::integrate_rotw(size_t i, size_t j) const
{
//#if defined __GNUC__
//#warning conjugate
//#endif
    complex<double> result = 0.0;
    for(size_t k = 0; k < tet_integration::gauss_num; k++)
        result += tet_integration::gauss_weights[k] *
                  rotw_pml(i, gauss_points_pml[k], gauss_points[k]).cj() *
                  rotw_pml(j, gauss_points_pml[k], gauss_points[k]).cj();
    return result * jacobian_pml;
}

complex<double> tetrahedron_pml::integrate_fw(cvector3(*func)(const point &, const phys_area &), size_t i) const
{
//#if defined __GNUC__
//#warning conjugate
//#endif
    complex<double> result = 0.0;
    for(size_t k = 0; k < tet_integration::gauss_num; k++)
        result += tet_integration::gauss_weights[k] *
                  func(gauss_points[k], get_phys_area()) *
                  w_pml(i, gauss_points_pml[k]).cj();
    return result * jacobian_pml;
}

cmatrix12 tetrahedron_pml::G() const
{
    cmatrix12 matr;
    for(size_t i = 0; i < 12; i++)
        for(size_t j = 0; j < 12; j++)
            matr[i][j] = 0.0;

    for(size_t i = 0; i < 6; i++)
        for(size_t j = 0; j <= i; j++)
            matr[j][i] = matr[i][j] = integrate_rotw(i, j);
    return matr;
}

cmatrix12 tetrahedron_pml::M() const
{
    cmatrix12 matr;
    for(size_t i = 0; i < 12; i++)
        for(size_t j = 0; j <= i; j++)
            matr[j][i] = matr[i][j] = integrate_w(i, j);
    return matr;
}

carray12 tetrahedron_pml::rp(cvector3(*func)(const point &, const phys_area &)) const
{
    carray12 arr;
    for(size_t i = 0; i < 12; i++)
        arr[i] = integrate_fw(func, i);
    return arr;
}
