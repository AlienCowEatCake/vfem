#if !defined(VFEM_H_INCLUDED)
#define VFEM_H_INCLUDED

#include "../common/common.h"
#include "../common/config.h"
#include "../elements/edge.h"
#include "../elements/triangle.h"
#include "../elements/tetrahedron.h"
#include "../vfem/phys.h"
#include "../vfem/slae.h"

#define VFEM_USE_PML

#if defined(VFEM_USE_PML)
typedef tetrahedron_pml finite_element;
#else
typedef tetrahedron finite_element;
#endif

// Правая часть
cvector3 func_rp(const point & p, phys_area & phys, void * data);

// Функция неоднородных первых краевых условий
cvector3 func_b1(const point & p, phys_area & phys, void * data);

// Функция аналитического решения
cvector3 func_true(const point & p, phys_area & phys, void * data);

#if defined(VFEM_USE_PML)
// Коэффициент S для PML
cvector3 get_s(const point & p, const tetrahedron_pml * fe, const phys_pml_area * phys_pml);

// Проверка, PML или нет
bool is_pml(const point & p, const finite_element * fe, const phys_pml_area * phys_pml);
#endif

// Класс векторный МКЭ
class VFEM
{
public:
    VFEM();

    // Конфигурация
    config_type config;

    // Ввод данных
    bool input_phys(const string & phys_filename);
    bool input_mesh(const string & gmsh_filename);
    // Процедура для сборки СЛАУ
    void make_struct();
    void make_data();
    // Запуск решения СЛАУ
    void solve();
    // Вывод данных в 3D сетке
    bool output(const string & tecplot_filename);
    // Вывод данных в 2D сетке
    bool output_slice(const string & tecplot_filename, char slice_var, double slice_val,
                      char var1, double min_var1, double max_var1, size_t num_var_1,
                      char var2, double min_var2, double max_var2, size_t num_var_2);
    // Вывод данных по линии
    bool output_line(const string & tecplot_filename, char line_var1, double line_val1,
                     char line_var2, double line_val2, char var3, double min_var3,
                     double max_var3, size_t num_var);

    // Поиск конечного элемента по точке
    finite_element * get_fe(const point & p) const;

    // Решение в точке
    cvector3 solution(const point & p) const;
    cvector3 solution(const point & p, const finite_element * fe) const;
    // Ротор решения в точке
    cvector3 rotor(const point & p) const;
    cvector3 rotor(const point & p, const finite_element * fe) const;
    // Дивергенция решения в точке
    complex<double> divergence(const point & p) const;
    complex<double> divergence(const point & p, const finite_element * fe) const;

    // Конечные элементы (тетраэдры)
    vector<finite_element> fes;

    void calculate_diff();

    // Основная СЛАУ
    SLAE slae;
    // СЛАУ на ядре
    SLAE ker_slae;
    // СЛАУ по границе
    SLAE surf_slae;

    // Узлы
    vector<point> nodes;
    // Ребра
    set<edge> edges;
    // Грани
    set<face> faces;
    // Ребра с источниками
    vector<edge_src> edges_src;
    // Треугольники
    vector<triangle *> trs;
    // Треугольники простые
    // (не используется при config.boundary_enabled == true)
    vector<triangle_base> trs_base;
    // Треугольники полные, для первых краевых
    // (не используется при config.boundary_enabled == false)
    vector<triangle_full> trs_full;
    // Точечные источники
    vector<pair<point, cvector3> > pss;
    // Физические области
    map<phys_id, phys_area> phys;

    // Число степеней свободы
    size_t dof_num;
    // Число степеней свободы ядра
    size_t ker_dof_num;
    // Соответствие глобальных степеней свободы и по границе
    // (не используется при config.boundary_enabled == false)
    map<size_t, size_t> global_to_local;
    // Степени свободы с первыми краевыми
    // (не используется при config.boundary_enabled == true)
    set<size_t> dof_first;
    // Степени свободы с первыми краевыми у ядра
    set<size_t> ker_dof_first;
    // Восьмиричное дерево поиска
    octree<finite_element> tree;

    // Получение количества степеней свободы тетраэдра
    size_t get_tet_dof_num(const tetrahedron_base * fe) const;
    // Получение степеней свободы тетраэдра в глобальной матрице
    size_t get_tet_dof(const tetrahedron_base * fe, size_t i) const;
    // Получение количества степеней свободы ядра тетраэдра
    size_t get_tet_ker_dof_num(const tetrahedron_base * fe) const;
    // Получение степеней свободы тетраэдра в матрице ядра
    size_t get_tet_ker_dof(const tetrahedron_base * fe, size_t i) const;
    // Получение степеней свободы треугольника в глобальной матрице
    size_t get_tr_dof(const triangle * tr, size_t i) const;
    // Получение степеней свободы треугольника в матрице ядра
    size_t get_tr_ker_dof(const triangle * tr, size_t i) const;
    // Получение степеней свободы треугольника в матрице по границе
    size_t get_tr_surf_dof(const triangle * tr, size_t i) const;

    // Генерация портрета глобальной матрицы
    void generate_portrait();
    // Генерация портрета глобальной матрицы ядра
    void generate_ker_portrait();
    // Сборка глобальной матрицы
    void assemble_matrix();
    // Применение краевых условий
    void applying_bound();
    // Генерация портрета по границе
    void generate_surf_portrait();
    // Применение источников на ребрах
    void apply_edges_sources();
    // Применение точечных источников
    void apply_point_sources();
    // Применение электродов
    void apply_electrodes();

protected:
    // Добавление локальной матрицы от одного КЭ в глобальную
    template<typename T>
    void process_fe_MpG(T * curr_fe);
    // Добавление локальной матрицы ядра от одного КЭ в глобальную матриицу ядра
    template<typename T>
    void process_fe_K(T * curr_fe);
    // Добавление локальной правой части от одного КЭ в глобальную правую часть
    template<typename T>
    void process_fe_rp(T * curr_fe);

    // Генерация абстрактного портрета
    template<typename U, typename V>
    void generate_abstract_portrait(size_t all_dof_num, const vector<U> & elems, SLAE & slae_curr,
                                    size_t (VFEM::*get_dof_num)(const V *) const,
                                    size_t (VFEM::*get_dof)(const V *, size_t) const) const;

    // Добавление ребра в множество ребер
    size_t add_edge(edge ed, set<edge> & edges_set);
    // Добавление грани в множество граней
    size_t add_face(face fc, set<face> & faces_set);

#if defined(VFEM_USE_PML)
    cpoint convert_point_to_pml(const point * p, const finite_element * fefe) const;
    void input_pml();
    // Параметры PML
    phys_pml_area phys_pml;
#endif

    // Проектирование на пространство ядра
    void to_kernel_space(const complex<double> * in, complex<double> * out) const;
    // Интерполяция на полное пространство
    void to_full_space(const complex<double> * in, complex<double> * out) const;
    // Скалярное произведение
    double dot_prod_self(const complex<double> * a) const;
    // Умножение матрицы с полного пространства на вектор
    void mul_matrix(const complex<double> * f, complex<double> * x) const;
    // Подсчет невязки
    void calc_residual(const complex<double> * x0, complex<double> * p) const;
};

#endif // VFEM_H_INCLUDED
