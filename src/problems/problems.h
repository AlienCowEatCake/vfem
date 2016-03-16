#if !defined PROBLEMS_H_INCLUDED
#define PROBLEMS_H_INCLUDED

#include "../common/common.h"
#include "../vfem/vfem.h"

// Контейнер области для сравнения
class diff_area
{
public:
    bool included;
    point p1, p2;
    diff_area()
    {
        included = true;
        p1 = point(-DBL_MAX, -DBL_MAX, -DBL_MAX);
        p2 = point(DBL_MAX, DBL_MAX, DBL_MAX);
    }
};

// Постпроцессор
void postprocessing(VFEM & v, const char * timebuf);

// Простое сравнение (на геометрически идентичных подобластях сеток и одинаковых базисах)
void compare_simple(VFEM & master, VFEM & slave, const vector<diff_area> & areas);

// Сложное сравнение (на произвольных подобластях сеток и базисах)
void compare_complex(VFEM & master, VFEM & slave, const vector<diff_area> & areas);

// Чтение конфигурационного файла для сравнения
bool input_diff(const string & filename, vector<diff_area> & areas);

#endif // PROBLEMS_H_INCLUDED
