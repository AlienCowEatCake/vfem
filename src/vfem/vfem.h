#if !defined VFEM_H_INCLUDED
#define VFEM_H_INCLUDED

#include "../common/common.h"
#include "../common/matrix.h"
#include "../geometry/point.h"
#include "../geometry/vector3.h"
#include "../elements/octal_tree.h"
#include "../elements/basis_config.h"
#include "../elements/edge.h"
#include "../elements/face.h"
#include "../elements/triangle.h"
#include "../elements/tetrahedron.h"
#include "../vfem/phys.h"
#include "../vfem/slae.h"

//#define VFEM_USE_PML
#define VFEM_USE_NONHOMOGENEOUS_FIRST
#define VFEM_USE_ANALYTICAL

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
    VFEM();
    ~VFEM();

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

    finite_element * get_fe(const point & p) const; // Поиск конечного элемента по точке

    cvector3 solution(const point & p) const;   // Решение в точке
    cvector3 rotor(const point & p) const;      // Ротор решения в точке
    cvector3 solution(const point & p, const finite_element * fe) const;
    cvector3 rotor(const point & p, const finite_element * fe) const;

    size_t fes_num;         // Число конечных элементов (тетраэдров)
    finite_element * fes;   // Конечные элементы (тетраэдры)

#if defined VFEM_USE_ANALYTICAL
    void calculate_diff() const;
#endif

    SLAE slae;              // Основная СЛАУ
#if defined VFEM_USE_NONHOMOGENEOUS_FIRST
    SLAE surf_slae;         // СЛАУ по границе
#endif
protected:
    size_t add_edge(edge ed, set<edge> & edges_set);    // Добавление ребра в множество ребер

    vector<point> nodes;    // Узлы

    size_t edges_num;   // Число ребер
    edge * edges;       // Ребра

#if BASIS_ORDER >= 2
    size_t faces_num;   // Число граней
    face * faces;       // Грани
#endif

    size_t edges_src_num;   // Число ребер с источниками
    edge_src * edges_src;   // Ребра с источниками

    map<phys_id, phys_area> phys;   // Физические области

    size_t trs_num;     // Число треугольников
    triangle * trs;     // Треугольники

#if defined VFEM_USE_NONHOMOGENEOUS_FIRST
    size_t dof_surf_num;        // Число степеней свободы с первыми неоднородными краевыми
    map<size_t, size_t> global_to_local;    // Соответствие глобальных степеней свободы и по границе
#else
    set<size_t> dof_first;      // Степени свободы с первыми краевыми
#endif
    size_t bound1_num;          // Число треугольников с первыми краевыми
    size_t bound2_num;          // Число треугольников со вторыми краевыми

    octal_tree<finite_element> tree;    // Восьмиричное дерево поиска

    void generate_portrait();       // Генерация портрета глобальной матрицы
    void assemble_matrix();         // Сборка глобальной матрицы
    void applying_bound();          // Применение краевых условий
#if defined VFEM_USE_NONHOMOGENEOUS_FIRST
    void generate_surf_portrait();  // Генерация портрета по границе
#endif
    void apply_edges_sources();     // Применение источников на ребрах
    void apply_point_sources();     // Применение точечных источников
    size_t pss_num;                 // Число точечных источников
    pair<point, cvector3> * pss;    // Точечные источники

#if defined VFEM_USE_PML
    cpoint convert_point_to_pml(const point * p, const finite_element * fefe) const;
    void input_pml();
    phys_pml_area phys_pml;         // Параметры PML
#endif
};

#endif // VFEM_H_INCLUDED
