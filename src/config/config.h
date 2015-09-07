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

// Вычислитель для вычисляемых значений
class evaluator
{
public:
    evaluator();
    cvector3 eval(const point & p, const phys_area * phys);
    map<size_t, array_t<parser<complex<double> >, 3> > values;
    array_t<parser<complex<double> >, 3> default_value;
};

// Класс конфигурации
class config_type
{
public:
    config_type();
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

    // ===== Boundary =====

    evaluator boundary;

    // ===== Right =====

    evaluator right;

    // ===== Analytical =====

    evaluator analytical;
    bool analytical_enabled;

protected:
    void load_defaults();
};

#endif // CONFIG_H

