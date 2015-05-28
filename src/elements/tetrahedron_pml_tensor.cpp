#include "tetrahedron.h"

void tetrahedron_pml_tensor::init_pml(cvector3(* get_s)(const point &, const tetrahedron_pml_tensor *, const phys_pml_area *))
{
    this->get_s = get_s;
}

complex<double> tetrahedron_pml_tensor::integrate_w(size_t i, size_t j) const
{
    using namespace tet_integration;
    matrix_t<complex<double>, 3, 3> lambda;
    for(size_t ik = 0; ik < 3; ik++)
        for(size_t jk = 0; jk < 3; jk++)
            lambda[ik][jk] = 0.0;
    complex<double> result = 0.0;
    for(size_t k = 0; k < gauss_num; k++)
    {
        cvector3 s = get_s(gauss_points[k], this, NULL);
        lambda[0][0] = s.y * s.z / s.x;
        lambda[1][1] = s.x * s.z / s.y;
        lambda[2][2] = s.x * s.y / s.z;
        result += gauss_weights[k] *
                  (lambda * cvector3(w(i, gauss_points[k]))) *
                  w(j, gauss_points[k]);
    }
    return result * jacobian;
}

complex<double> tetrahedron_pml_tensor::integrate_rotw(size_t i, size_t j) const
{
    using namespace tet_integration;
    matrix_t<complex<double>, 3, 3> lambda1;
    for(size_t ik = 0; ik < 3; ik++)
        for(size_t jk = 0; jk < 3; jk++)
            lambda1[ik][jk] = 0.0;
    complex<double> result = 0.0;
    for(size_t k = 0; k < gauss_num; k++)
    {
        cvector3 s = get_s(gauss_points[k], this, NULL);
        lambda1[0][0] = s.x / (s.y * s.z);
        lambda1[1][1] = s.y / (s.x * s.z);
        lambda1[2][2] = s.z / (s.x * s.y);
        result += gauss_weights[k] *
                  (lambda1 * cvector3(rotw(i, gauss_points[k]))) *
                  rotw(j, gauss_points[k]);
    }
    return result * jacobian;
}

complex<double> tetrahedron_pml_tensor::integrate_fw(cvector3(*func)(const point &, const phys_area &), size_t i) const
{
    using namespace tet_integration;
    complex<double> result = 0.0;
    for(size_t k = 0; k < gauss_num; k++)
        result += gauss_weights[k] *
                  func(gauss_points[k], get_phys_area()) *
                  w(i, gauss_points[k]);
    return result * jacobian;
}

matrix_t<complex<double>, basis::tet_bf_num, basis::tet_bf_num>
tetrahedron_pml_tensor::G() const
{
    using namespace basis;
    matrix_t<complex<double>, tet_bf_num, tet_bf_num> matr;
    for(size_t i = 0; i < tet_bf_num; i++)
        for(size_t j = 0; j <= i; j++)
            matr[j][i] = matr[i][j] = integrate_rotw(i, j);
    return matr;
}

matrix_t<complex<double>, basis::tet_bf_num, basis::tet_bf_num>
tetrahedron_pml_tensor::M() const
{
    using namespace basis;
    matrix_t<complex<double>, tet_bf_num, tet_bf_num> matr;
    for(size_t i = 0; i < tet_bf_num; i++)
        for(size_t j = 0; j <= i; j++)
            matr[j][i] = matr[i][j] = integrate_w(i, j);
    return matr;
}

array_t<complex<double>, basis::tet_bf_num>
tetrahedron_pml_tensor::rp(cvector3(*func)(const point &, const phys_area &)) const
{
    using namespace basis;
    array_t<complex<double>, tet_bf_num> arr;
    for(size_t i = 0; i < tet_bf_num; i++)
        arr[i] = integrate_fw(func, i);
    return arr;
}
