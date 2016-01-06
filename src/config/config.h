#ifndef CONFIG_H
#define CONFIG_H

#include <cstdlib>
#include <string>
#include <map>
#include "../config/parser.h"
#include "../vfem/phys.h"
#include "../geometry/point.h"
#include "../geometry/vector3.h"
#include "../common/matrix.h"

using namespace std;

//                                  BASIS_ORDER BASIS_TYPE
// Первый порядок I типа (неполный)      1           1
// Первый порядок II типа (полный)       1           2
// Второй порядок I типа (неполный)      2           1
// Второй порядок II типа (полный)       2           2

// Конфигурация базиса
struct basis_type
{
    // Количество базисных функций на тетраэдрах
    size_t tet_bf_num;
    // Количество базисных функций ядра на тетраэдрах
    size_t tet_ker_bf_num;
    // Количество базисных функций на треугольниках
    size_t tr_bf_num;
    // Количество базисных функций ядра на треугольниках
    size_t tr_ker_bf_num;
    // Порядок базиса
    size_t order;
    size_t type;
};

// Вычислитель для вычисляемых значений
class evaluator
{
public:
    evaluator();
    // Вычислители по физическим областям
    map<size_t, array_t<parser<complex<double> >, 3> > values;
    // Вычислители по умолчанию
    array_t<parser<complex<double> >, 3> default_value;
    // Тип JIT-компилятора в вычислителях
    enum jit_types
    {
        JIT_DISABLE,
        JIT_INLINE,
        JIT_EXTCALL
    };
};

// Конфигурация 1D постпроцессора
struct postprocessor_1d
{
    char line_var1_name;
    double line_var1_value;
    char line_var2_name;
    double line_var2_value;
    char var_name;
    double var_from;
    double var_to;
    size_t var_num;
};

// Конфигурация 2D постпроцессора
struct postprocessor_2d
{
    char slice_var_name;
    double slice_var_value;
    char var1_name;
    double var1_from;
    double var1_to;
    size_t var1_num;
    char var2_name;
    double var2_from;
    double var2_to;
    size_t var2_num;
};

// Общая конфигурация постпроцессора
class postprocessor
{
public:
    unsigned char type;
    string filename;
    bool timestamp;
    union
    {
        postprocessor_1d param_1d;
        postprocessor_2d param_2d;
    };
    postprocessor() : type(3), filename("plot.plt"), timestamp(false)
    {
        memset(&param_1d, 0, sizeof(postprocessor_1d));
        memset(&param_2d, 0, sizeof(postprocessor_2d));
    }
};

// Класс конфигурации
class config_type
{
public:
    config_type();
    // Загрузка значений из файла
    bool load(const string & filename);

    // ===== VFEM =====

    // Конфигурация базиса
    basis_type basis;
    // Точность решения СЛАУ
    double eps_slae;
    // Точность решения СЛАУ для краевых
    double eps_slae_bound;
    // Точность начального решения перед запуском v-цикла
    double gamma_v_cycle_0;
    // Точность решения на полном пространстве
    double gamma_v_cycle_full;
    // Точность решения на пространстве ядра
    double gamma_v_cycle_ker;
    // Максимальное локальное число итераций v-цикла
    size_t max_iter_v_cycle_local;
    // Сетка
    string filename_mesh;
    // Параметры физических областей
    string filename_phys;
    // Куда сохранять веса решения
    string filename_slae;
    // Тип JIT-компилятора в вычислителях
    evaluator::jit_types jit_type;

    // ===== Boundary =====

    evaluator boundary;
    bool boundary_enabled;

    // ===== Right =====

    evaluator right;
    bool right_enabled;

    // ===== Analytical =====

    evaluator analytical;
    bool analytical_enabled;

    // ===== Postprocessing =====

    map<size_t, postprocessor> post;

protected:
    // Загрузка значений по-умолчанию
    void load_defaults();
    // Пост-загрузочная инициализация
    bool init(bool status);
};

#endif // CONFIG_H

