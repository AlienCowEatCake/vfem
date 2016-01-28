#if !defined TETRAHEDRON_H_INCLUDED
#define TETRAHEDRON_H_INCLUDED

#include "../common/common.h"
#include "../common/matrix.h"
#include "../common/cubatures.h"
#include "../common/trio.h"
#include "../config/config.h"
#include "../geometry/point.h"
#include "../geometry/vector3.h"
#include "../elements/edge.h"
#include "../elements/face.h"
#include "../vfem/phys.h"

using namespace tet_integration_8;

typedef cvector3(* eval_func)(const point &, const phys_area &, void *);

// Индексы для построения базисных функций на тетраэдрах
namespace tet_basis_indexes
{
    // Edges (Ребра) // k, l : k < l
    extern const size_t ind_e[6][2];
    // Faces (Грани) // j, k, l : j < k < l
    extern const size_t ind_f[4][3];
}

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
    face * faces[4];    // Грани
    const face & get_face(size_t i) const;

    phys_area * phys;   // Физическая область
    const phys_area & get_phys_area() const;

    const basis_type * basis;   // Параметры базиса

    // Базисные функции
    vector3 w(size_t i, const point & p) const;
    // Роторы базисных функций
    vector3 rotw(size_t i, const point & p) const;
    // Базисные функции ядра
    vector3 kerw(size_t i, const point & p) const;

    point barycenter;

    bool in_cube(double x0, double x1, double y0, double y1, double z0, double z1) const;

    trio<double, vector3, cvector3>
    diff_normL2(const array_t<complex<double> > & q, eval_func func, void * data) const;
    trio<double, vector3, cvector3>
    diff_normL2(const array_t<complex<double> > & q, const array_t<complex<double> > & q_true) const;
    trio<double, vector3, cvector3>
    normL2(eval_func func, void * data) const;
    trio<double, vector3, cvector3>
    normL2(const array_t<complex<double> > & q_true) const;

protected:
    // Матрица L-координат
    matrix_t<double, 4, 4> L;
    // Градиент L-координаты
    vector3 grad_lambda(size_t i) const;
    // L-координаты
    double lambda(size_t i, const point & p) const;

    // Точки Гаусса
    point gauss_points[tet_integration::gauss_num];
    // Якобиан
    double jacobian;

    // Параметры прямых для дерева
    double edges_a[6][3], edges_b[6][3];
};

// Класс тетраэдр (обычный)
class tetrahedron : public tetrahedron_base
{
public:
    // Локальная матрица жескости
    matrix_t<double> G() const;
    // Локальная матрица массы
    matrix_t<double> M() const;
    // Локальная правая часть
    array_t<complex<double> > rp(eval_func func, void * data) const;
    // Локальная матрица ядра
    matrix_t<double> K() const;

protected:
    // Интеграл от бф
    double integrate_w(size_t i, size_t j) const;
    // Интеграл от ротора бф
    double integrate_rotw(size_t i, size_t j) const;
    complex<double> integrate_fw(eval_func func, size_t i, void * data) const;
    // Интегралы от базисных функций ядра
    double integrate_kerw(size_t i, size_t j) const;
};

// Класс тетраэдр (для работы с PML-краевыми)
class tetrahedron_pml : public tetrahedron_base
{
public:
    void init_pml(cvector3(* get_s)(const point &, const tetrahedron_pml *, const phys_pml_area *), const phys_pml_area * phys_pml, const cpoint * nodes_pml);

    // Локальная матрица жескости
    matrix_t<complex<double> > G() const;
    // Локальная матрица массы
    matrix_t<complex<double> > M() const;
    // Локальная правая часть
    array_t<complex<double> > rp(eval_func func, void * data) const;
    // Локальная матрица ядра
    matrix_t<complex<double> > K() const;

protected:
    cvector3(* get_s)(const point &, const tetrahedron_pml *, const phys_pml_area *);
    const phys_pml_area * phys_pml;

    // Матрица L-координат (в PML)
    matrix_t<complex<double>, 4, 4> L_pml;
    // Градиент L-координаты (в PML)
    cvector3 grad_lambda_pml(size_t i) const;
    // L-координаты (в PML)
    complex<double> lambda_pml(size_t i, const cpoint & p) const;

    // Точки Гаусса (в PML)
    cpoint gauss_points_pml[tet_integration::gauss_num];
    // Якобиан (в PML)
    complex<double> jacobian_pml;

    // Базисные функции (в PML)
    cvector3 w_pml(size_t i, const cpoint & p) const;
    // Роторы базисных функций (в PML)
    cvector3 rotw_pml(size_t i, const cpoint & p, const point & p_non_PML) const;
    // Базисные функции ядра (в PML)
    cvector3 kerw_pml(size_t i, const cpoint & p, const point & p_non_PML) const;

    // Интеграл от бф
    complex<double> integrate_w(size_t i, size_t j) const;
    // Интеграл от ротора бф
    complex<double> integrate_rotw(size_t i, size_t j) const;
    complex<double> integrate_fw(eval_func func, size_t i, void * data) const;
    // Интегралы от базисных функций ядра
    complex<double> integrate_kerw(size_t i, size_t j) const;
};

#endif // TETRAHEDRON_H_INCLUDED
