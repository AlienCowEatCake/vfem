#ifndef BASIS_H_INCLUDED
#define BASIS_H_INCLUDED

#include <cstring>
#include <cmath>

//                                  BASIS_ORDER BASIS_TYPE
// Первый порядок I типа (неполный)      1           1
// Первый порядок II типа (полный)       1           2
// Второй порядок I типа (неполный)      2           1
// Второй порядок II типа (полный)       2           2
#define BASIS_ORDER 1
#define BASIS_TYPE  2

// ====== Базис первого неполного порядка =====================================

#if BASIS_ORDER == 1 && BASIS_TYPE == 1
namespace basis
{
    static const size_t tet_bf_num = 6;
    static const size_t tr_bf_num = 3;
}
#endif

// ====== Базис первого полного порядка =======================================

#if BASIS_ORDER == 1 && BASIS_TYPE == 2
namespace basis
{
    static const size_t tet_bf_num = 12;
    static const size_t tr_bf_num = 6;
}
#endif

// ====== Базис второго неполного порядка =====================================

#if BASIS_ORDER == 2 && BASIS_TYPE == 1
namespace basis
{
    static const size_t tet_bf_num = 20;
    static const size_t tr_bf_num = 8;
}
#endif

// ====== Базис второго полного порядка =======================================

#if BASIS_ORDER == 2 && BASIS_TYPE == 2
namespace basis
{
    static const size_t tet_bf_num = 30;
    static const size_t tr_bf_num = 12;
}
#endif

// ============================================================================

// Индексы для построения базисных функций на тетраэдрах
namespace tet_basis_indexes
{
    // Edges (Ребра) // k, l : k < l
    extern const size_t ind_e[6][2];
    // Faces (Грани) // j, k, l : j < k < l
    extern const size_t ind_f[4][3];
}

// ============================================================================

// Индексы для построения базисных функций на треугольниках
namespace tr_basis_indexes
{
    // Edges (Ребра) // k, l : k < l
    extern const size_t ind_e[3][2];
}

// ============================================================================

// Численное интегрирование на тетраэдрах
// https://people.fh-landshut.de/~maurer/numeth/node148.html

namespace tet_integration_2
{
    namespace tet_integration
    {
        static const size_t gauss_num = 4;
        extern const double gauss_weights[gauss_num];
        extern const double gauss_points_master[gauss_num][4];
    }
}

namespace tet_integration_3
{
    namespace tet_integration
    {
        static const size_t gauss_num = 5;
        extern const double gauss_weights[gauss_num];
        extern const double gauss_points_master[gauss_num][4];
    }
}

namespace tet_integration_4
{
    namespace tet_integration
    {
        static const size_t gauss_num = 11;
        extern const double gauss_weights[gauss_num];
        extern const double gauss_points_master[gauss_num][4];
    }
}

// ============================================================================

// Численное интегрирование на треугольниках
// https://people.fh-landshut.de/~maurer/numeth/node147.html

namespace tr_integration_2
{
    namespace tr_integration
    {
        static const size_t gauss_num = 3;
        extern const double gauss_weights[gauss_num];
        extern const double gauss_points_master[gauss_num][3];
    }
}

namespace tr_integration_3
{
    namespace tr_integration
    {
        static const size_t gauss_num = 4;
        extern const double gauss_weights[gauss_num];
        extern const double gauss_points_master[gauss_num][3];
    }
}

namespace tr_integration_3_2
{
    namespace tr_integration
    {
        static const size_t gauss_num = 7;
        extern const double gauss_weights[gauss_num];
        extern const double gauss_points_master[gauss_num][3];
    }
}

namespace tr_integration_5
{
    namespace tr_integration
    {
        static const size_t gauss_num = 7;
        extern const double gauss_weights[gauss_num];
        extern const double gauss_points_master[gauss_num][3];
    }
}

// ============================================================================

// Интегрирование особо высоких порядков
// http://lsec.cc.ac.cn/~tcui/myinfo/paper/quad.pdf

namespace tet_integration_14
{
    namespace tet_integration
    {
        static const size_t gauss_num = 236;
        extern const double gauss_weights[gauss_num];
        extern const double gauss_points_master[gauss_num][4];
    }
}

namespace tr_integration_21
{
    namespace tr_integration
    {
        static const size_t gauss_num = 91;
        extern const double gauss_weights[gauss_num];
        extern const double gauss_points_master[gauss_num][3];
    }
}

// ============================================================================

#if BASIS_ORDER == 1
using namespace tet_integration_4;
using namespace tr_integration_5;
#endif

#if BASIS_ORDER == 2
using namespace tet_integration_14;
using namespace tr_integration_21;
#endif

#endif // BASIS_H_INCLUDED
