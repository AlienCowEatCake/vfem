#include "basis_config.h"

// ============================================================================

// Индексы для построения базисных функций на тетраэдрах
namespace tet_basis_indexes
{
    // Edges (Ребра) // k, l : k < l
    const size_t ind_e[6][2] =
    {
        { 0, 1 },
        { 0, 2 },
        { 0, 3 },
        { 1, 2 },
        { 1, 3 },
        { 2, 3 }
    };
    // Faces (Грани) // j, k, l : j < k < l
    const size_t ind_f[4][3] =
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
    const size_t ind_e[3][2] =
    {
        { 0, 1 },
        { 0, 2 },
        { 1, 2 }
    };
}

// ============================================================================
