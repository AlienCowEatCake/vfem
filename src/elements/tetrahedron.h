#if !defined TETRAHEDRON_H_INCLUDED
#define TETRAHEDRON_H_INCLUDED

#include "../common/common.h"
#include "../common/matrix.h"
#include "../geometry/point.h"
#include "../geometry/vector3.h"
#include "../elements/edge.h"
#include "../vfem/phys.h"
#include "../elements/basis_config.h"

// Класс тетраэдр (абстрактный)
class tetrahedron_base
{
public:
    tetrahedron_base();
    void init();
    bool inside(const point & p) const;
    bool inside(double x, double y, double z) const;

    point * nodes[4];   // Узлы
    edge * edges[6];    // Ребра
    const point & get_node(size_t i) const;
    const edge & get_edge(size_t i) const;

    phys_area * phys;    // Физическая область
    const phys_area & get_phys_area() const;    // Получение физической области

    vector3 w(size_t i, const point & p) const;     // Базисные функции
    vector3 rotw(size_t i, const point & p) const;  // Роторы базисных функций

    point barycenter;

    bool inside_tree(double x0, double x1, double y0, double y1, double z0, double z1) const;
    double edges_a[6][3], edges_b[6][3];

    double diff_normL2(const carray12 & q, cvector3(*func)(const point &)) const;
    double normL2(cvector3(*func)(const point &)) const;

protected:
    matrix_t<double, 4, 4> L;   // Матрица L-координат
    vector3 grad_lambda(size_t i) const;    // градиент L-координаты
    double lambda(size_t i, const point & p) const; // L-координаты

    point gauss_points[tet_integration::gauss_num]; // Точки Гаусса
    double jacobian;            // Якобиан
};

// Класс тетраэдр (обычный)
class tetrahedron : public tetrahedron_base
{
public:
    // Локальная матрица жескости
    matrix_t<double, basis::tet_bf_num, basis::tet_bf_num> G() const;
    // Локальная матрица массы
    matrix_t<double, basis::tet_bf_num, basis::tet_bf_num> M() const;
    // Локальная правая часть
    array_t<complex<double>, basis::tet_bf_num> rp(cvector3(*func)(const point & , const phys_area &)) const;

protected:
    double integrate_w(size_t i, size_t j) const;       // Интеграл от бф
    double integrate_rotw(size_t i, size_t j) const;    // Интеграл от ротора бф
    complex<double> integrate_fw(cvector3(*func)(const point & , const phys_area &), size_t i) const;
};

// Класс тетраэдр (для работы с PML-краевыми)
class tetrahedron_pml : public tetrahedron_base
{
public:
    void init_pml(cvector3(* get_s)(const point &, const tetrahedron_pml *, const phys_pml_area *), const phys_pml_area * phys_pml);

    // Локальная матрица жескости
    matrix_t<complex<double>, basis::tet_bf_num, basis::tet_bf_num> G() const;
    // Локальная матрица массы
    matrix_t<complex<double>, basis::tet_bf_num, basis::tet_bf_num> M() const;
    // Локальная правая часть
    array_t<complex<double>, basis::tet_bf_num> rp(cvector3(*func)(const point & , const phys_area &)) const;

    cpoint * nodes_pml[4];   // Узлы (в PML)
    const cpoint & get_node_pml(size_t i) const;

protected:
    cvector3(* get_s)(const point &, const tetrahedron_pml *, const phys_pml_area *);
    const phys_pml_area * phys_pml;

    matrix_t<complex<double>, 4, 4> L_pml;  // Матрица L-координат (в PML)
    cvector3 grad_lambda_pml(size_t i) const;    // градиент L-координаты (в PML)
    complex<double> lambda_pml(size_t i, const cpoint & p) const; // L-координаты (в PML)

    cpoint gauss_points_pml[tet_integration::gauss_num];    // Точки Гаусса (в PML)
    complex<double> jacobian_pml;   // Якобиан (в PML)

    cvector3 w_pml(size_t i, const cpoint & p) const;     // Базисные функции (в PML)
    cvector3 rotw_pml(size_t i, const cpoint & p, const point & p_non_PML) const;  // Роторы базисных функций (в PML)

    complex<double> integrate_w(size_t i, size_t j) const;       // Интеграл от бф
    complex<double> integrate_rotw(size_t i, size_t j) const;    // Интеграл от ротора бф
    complex<double> integrate_fw(cvector3(*func)(const point & , const phys_area &), size_t i) const;
};

#endif // TETRAHEDRON_H_INCLUDED
