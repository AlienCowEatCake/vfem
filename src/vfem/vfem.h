#if !defined VFEM_H_INCLUDED
#define VFEM_H_INCLUDED

#include "../common/common.h"
#include "../common/matrix.h"
#include "../common/basis_config.h"
#include "../geometry/point.h"
#include "../geometry/vector3.h"
#include "../elements/octal_tree.h"
#include "../elements/edge.h"
#include "../elements/face.h"
#include "../elements/triangle.h"
#include "../elements/tetrahedron.h"
#include "../vfem/phys.h"
#include "../vfem/slae.h"

#define VFEM_USE_PML
//#define VFEM_USE_NONHOMOGENEOUS_FIRST
//#define VFEM_USE_ANALYTICAL

#if defined VFEM_USE_PML
typedef tetrahedron_pml finite_element;
typedef matrix_t<complex<double>, basis::tet_bf_num, basis::tet_bf_num> l_matrix;
#else
typedef tetrahedron finite_element;
typedef matrix_t<double, basis::tet_bf_num, basis::tet_bf_num> l_matrix;
#endif

#if defined VFEM_USE_NONHOMOGENEOUS_FIRST
typedef triangle_full triangle;
#else
typedef triangle_base triangle;
#endif

// Правая часть
cvector3 func_rp(const point & p, const phys_area & phys);

// Функция неоднородных первых краевых условий
#if defined VFEM_USE_NONHOMOGENEOUS_FIRST
cvector3 func_b1(const point & p, const triangle * tr);
#endif

// Функция аналитического решения
#if defined VFEM_USE_ANALYTICAL
cvector3 func_true(const point & p);
#endif

// Коэффициент S для PML
#if defined VFEM_USE_PML
cvector3 get_s(const point & p, const tetrahedron_pml * fe, const phys_pml_area * phys_pml);
#endif

// Проверка, PML или нет
#if defined VFEM_USE_PML
bool is_pml(const point & p, const finite_element * fe);
#endif

// Класс векторный МКЭ
class VFEM
{
public:
    // Ввод данных
    void input_phys(const string & phys_filename);
    void input_mesh(const string & gmsh_filename);
    // Процедура решения
    void solve();
    // Вывод данных в 3D сетке
    void output(const string & tecplot_filename);
    // Вывод данных в 2D сетке
    void output_slice(const string & tecplot_filename, char slice_var, double slice_val,
                      char var1, double min_var1, double max_var1, double step_var1,
                      char var2, double min_var2, double max_var2, double step_var2);

    // Поиск конечного элемента по точке
    finite_element * get_fe(const point & p) const;

    // Решение в точке
    cvector3 solution(const point & p) const;
    cvector3 solution(const point & p, const finite_element * fe) const;
    // Ротор решения в точке
    cvector3 rotor(const point & p) const;
    cvector3 rotor(const point & p, const finite_element * fe) const;

    // Конечные элементы (тетраэдры)
    vector<finite_element> fes;

#if defined VFEM_USE_ANALYTICAL
    void calculate_diff() const;
#endif

    // Основная СЛАУ
    SLAE slae;
#if defined VFEM_USE_NONHOMOGENEOUS_FIRST
    // СЛАУ по границе
    SLAE surf_slae;
#endif
protected:
    // Добавление ребра в множество ребер
    size_t add_edge(edge ed, set<edge> & edges_set);
#if BASIS_ORDER >= 2
    // Добавление грани в множество граней
    size_t add_face(face fc, set<face> & faces_set);
#endif

    // Узлы
    vector<point> nodes;
    // Ребра
    set<edge> edges;
#if BASIS_ORDER >= 2
    // Грани
    set<face> faces;
#endif
    // Ребра с источниками
    vector<edge_src> edges_src;
    // Треугольники
    vector<triangle> trs;
    // Точечные источники
    vector<pair<point, cvector3> > pss;
    // Физические области
    map<phys_id, phys_area> phys;

    // Число степеней свободы
    size_t dof_num;
#if defined VFEM_USE_NONHOMOGENEOUS_FIRST
    // Соответствие глобальных степеней свободы и по границе
    map<size_t, size_t> global_to_local;
#else
    // Степени свободы с первыми краевыми
    set<size_t> dof_first;
#endif
    // Восьмиричное дерево поиска
    octal_tree<finite_element> tree;

    // Генерация портрета глобальной матрицы
    void generate_portrait();
    // Сборка глобальной матрицы
    void assemble_matrix();
    // Применение краевых условий
    void applying_bound();
#if defined VFEM_USE_NONHOMOGENEOUS_FIRST
    // Генерация портрета по границе
    void generate_surf_portrait();
#endif
    // Применение источников на ребрах
    void apply_edges_sources();
    // Применение точечных источников
    void apply_point_sources();

#if defined VFEM_USE_PML
    cpoint convert_point_to_pml(const point * p, const finite_element * fefe) const;
    void input_pml();
    // Параметры PML
    phys_pml_area phys_pml;
#endif
};

#endif // VFEM_H_INCLUDED
