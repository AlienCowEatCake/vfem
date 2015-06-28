#ifndef BASIS_H_INCLUDED
#define BASIS_H_INCLUDED

#include <cstring>
#include <cmath>
#include "cubatures.h"

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
    static const size_t tet_ker_bf_num = 4;
    static const size_t tr_bf_num = 3;
    static const size_t tr_ker_bf_num = 3;
}
#endif

// ====== Базис первого полного порядка =======================================

#if BASIS_ORDER == 1 && BASIS_TYPE == 2
namespace basis
{
    static const size_t tet_bf_num = 12;
    static const size_t tet_ker_bf_num = 10;
    static const size_t tr_bf_num = 6;
    static const size_t tr_ker_bf_num = 6;
}
#endif

// ====== Базис второго неполного порядка =====================================

#if BASIS_ORDER == 2 && BASIS_TYPE == 1
namespace basis
{
    static const size_t tet_bf_num = 20;
    static const size_t tet_ker_bf_num = 10;
    static const size_t tr_bf_num = 8;
    static const size_t tr_ker_bf_num = 6;
}
#endif

// ====== Базис второго полного порядка =======================================

#if BASIS_ORDER == 2 && BASIS_TYPE == 2
namespace basis
{
    static const size_t tet_bf_num = 30;
    static const size_t tet_ker_bf_num = 20;
    static const size_t tr_bf_num = 12;
    static const size_t tr_ker_bf_num = 10;
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

#if BASIS_ORDER == 1
using namespace tet_integration_4;
using namespace tr_integration_5;
#endif

#if BASIS_ORDER == 2
using namespace tet_integration_8;
using namespace tr_integration_8;
#endif

#endif // BASIS_H_INCLUDED
