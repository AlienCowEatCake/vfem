#ifndef BASIS_H_INCLUDED
#define BASIS_H_INCLUDED

#include <cstring>
#include <cmath>
#include "integration.h"

//                                  BASIS_ORDER BASIS_TYPE
// Первый порядок I типа (неполный)      1           1
// Первый порядок II типа (полный)       1           2
// Второй порядок I типа (неполный)      2           1
// Второй порядок II типа (полный)       2           2
#define BASIS_ORDER 2
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
    static const size_t ind_e[6][2] =
    {
        { 0, 1 },
        { 0, 2 },
        { 0, 3 },
        { 1, 2 },
        { 1, 3 },
        { 2, 3 }
    };
    // Faces (Грани) // j, k, l : j < k < l
    static const size_t ind_f[4][3] =
    {
        { 0, 1, 2 },
        { 0, 1, 3 },
        { 0, 2, 3 },
        { 1, 2, 3 }
    };
}

// ============================================================================

// Индексы для построения базисных функций на треугольниках
namespace tr_basis_indexes
{
    // Edges (Ребра) // k, l : k < l
    static const size_t ind_e[3][2] =
    {
        { 0, 1 },
        { 0, 2 },
        { 1, 2 }
    };
}

// ============================================================================

#endif // BASIS_H_INCLUDED
