#ifndef CONFIG_H
#define CONFIG_H

#include <cstdlib>
#include <string>
#include <map>
#include "parser/parser.h"
#include "../vfem/phys.h"
#include "../geometry/point.h"
#include "../geometry/vector3.h"

using namespace std;

// Вычислитель для вычисляемых значений
class evaluator
{
public:
    evaluator();
    bool parse(const string & str);
    cvector3 eval(const point & p, const phys_area * phys);

protected:
    map<size_t, parser[3]> values;
    parser default_value[3];
};

// Класс конфигурации
class config
{
public:
    bool load(const string & filename);

    // ===== VFEM =====

    // Порядок базиса
    size_t basis_order;
    size_t basis_type;
    // Точность решения СЛАУ
    double eps_slae;
    // Точность решения СЛАУ для краевых
    double eps_slae_bound;
    // Точность решения на полном пространстве
    double gamma_v_cycle_full;
    // Точность решения на пространстве ядра
    double gamma_v_cycle_ker;
    // Сетка
    string filename_mesh;
    // Параметры физических областей
    string filename_phys;
    // Куда сохранять веса решения
    string filename_slae;

protected:
    void load_defaults();
};

#endif // CONFIG_H

