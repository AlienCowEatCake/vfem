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

// Численное интегрирование на тетраэдрах
// https://people.fh-landshut.de/~maurer/numeth/node148.html

namespace tet_integration_2
{
    namespace tet_integration
    {
        const double gauss_weights[gauss_num] =
        {
            1.0 / 24.0,
            1.0 / 24.0,
            1.0 / 24.0,
            1.0 / 24.0
        };
        const double gauss_a = (5.0 - sqrt(5.0)) / 20.0;
        const double gauss_b = (5.0 + 3.0 * sqrt(5.0)) / 20.0;
        const double gauss_points_master[gauss_num][4] =
        {
            { 1.0 - gauss_b - 2.0 * gauss_a, gauss_b, gauss_a, gauss_a },
            { 1.0 - gauss_b - 2.0 * gauss_a, gauss_a, gauss_b, gauss_a },
            { 1.0 - gauss_b - 2.0 * gauss_a, gauss_a, gauss_a, gauss_b },
            { 1.0 - 3.0 * gauss_a,           gauss_a, gauss_a, gauss_a }
        };
    }
}

namespace tet_integration_3
{
    namespace tet_integration
    {
        const double gauss_weights[gauss_num] =
        {
            - 2.0 / 15.0,
            3.0 / 40.0,
            3.0 / 40.0,
            3.0 / 40.0,
            3.0 / 40.0
        };
        const double gauss_a = 1.0 / 4.0;
        const double gauss_b = 1.0 / 2.0;
        const double gauss_c = 1.0 / 6.0;
        const double gauss_points_master[gauss_num][4] =
        {
            { 1.0 - 3.0 * gauss_a,           gauss_a, gauss_a, gauss_a },
            { 1.0 - gauss_b - 2.0 * gauss_c, gauss_b, gauss_c, gauss_c },
            { 1.0 - gauss_b - 2.0 * gauss_c, gauss_c, gauss_b, gauss_c },
            { 1.0 - gauss_b - 2.0 * gauss_c, gauss_c, gauss_c, gauss_b },
            { 1.0 - 3.0 * gauss_c,           gauss_c, gauss_c, gauss_c }
        };
    }
}

namespace tet_integration_4
{
    namespace tet_integration
    {
        const double gauss_weights[gauss_num] =
        {
            - 74.0 / 5625.0,
            343.0 / 45000.0,
            343.0 / 45000.0,
            343.0 / 45000.0,
            343.0 / 45000.0,
            56.0 / 2250.0,
            56.0 / 2250.0,
            56.0 / 2250.0,
            56.0 / 2250.0,
            56.0 / 2250.0,
            56.0 / 2250.0
        };
        const double gauss_a = 1.0 / 4.0;
        const double gauss_b = 11.0 / 14.0;
        const double gauss_c = 5.0 / 70.0;
        const double gauss_d = (1.0 + sqrt(5.0 / 14.0)) / 4.0;
        const double gauss_e = (1.0 - sqrt(5.0 / 14.0)) / 4.0;
        const double gauss_points_master[gauss_num][4] =
        {
            { 1.0 - 3.0 * gauss_a,           gauss_a, gauss_a, gauss_a },
            { 1.0 - gauss_b - 2.0 * gauss_c, gauss_b, gauss_c, gauss_c },
            { 1.0 - gauss_b - 2.0 * gauss_c, gauss_c, gauss_b, gauss_c },
            { 1.0 - gauss_b - 2.0 * gauss_c, gauss_c, gauss_c, gauss_b },
            { 1.0 - 3.0 * gauss_c,           gauss_c, gauss_c, gauss_c },
            { 1.0 - gauss_d - 2.0 * gauss_e, gauss_d, gauss_e, gauss_e },
            { 1.0 - gauss_d - 2.0 * gauss_e, gauss_e, gauss_d, gauss_e },
            { 1.0 - gauss_d - 2.0 * gauss_e, gauss_e, gauss_e, gauss_d },
            { 1.0 - gauss_e - 2.0 * gauss_d, gauss_e, gauss_d, gauss_d },
            { 1.0 - gauss_e - 2.0 * gauss_d, gauss_d, gauss_e, gauss_d },
            { 1.0 - gauss_e - 2.0 * gauss_d, gauss_d, gauss_d, gauss_e }
        };
    }
}

// ============================================================================

// Численное интегрирование на треугольниках
// https://people.fh-landshut.de/~maurer/numeth/node147.html

namespace tr_integration_2
{
    namespace tr_integration
    {
        const double gauss_weights[gauss_num] =
        {
            1.0 / 6.0,
            1.0 / 6.0,
            1.0 / 6.0
        };
        const double gauss_points_master[gauss_num][3] =
        {
            { 1.0 / 2.0, 1.0 / 2.0, 0.0       },
            { 0.0,       1.0 / 2.0, 1.0 / 2.0 },
            { 1.0 / 2.0, 0.0,       1.0 / 2.0 }
        };
    }
}

namespace tr_integration_3
{
    namespace tr_integration
    {
        const double gauss_weights[gauss_num] =
        {
            - 9.0 / 32.0,
            25.0 / 96.0,
            25.0 / 96.0,
            25.0 / 96.0
        };
        const double gauss_a = 1.0 / 3.0;
        const double gauss_b = 3.0 / 5.0;
        const double gauss_c = 1.0 / 5.0;
        const double gauss_points_master[gauss_num][3] =
        {
            { 1.0 - 2.0 * gauss_a,     gauss_a, gauss_a },
            { 1.0 - gauss_b - gauss_c, gauss_b, gauss_c },
            { 1.0 - gauss_b - gauss_c, gauss_c, gauss_b },
            { 1.0 - 2.0 * gauss_c,     gauss_c, gauss_c }
        };
    }
}

namespace tr_integration_3_2
{
    namespace tr_integration
    {
        const double gauss_weights[gauss_num] =
        {
            9.0 / 40.0,
            1.0 / 15.0,
            1.0 / 15.0,
            1.0 / 15.0,
            1.0 / 40.0,
            1.0 / 40.0,
            1.0 / 40.0
        };
        const double gauss_a = 1.0 / 3.0;
        const double gauss_b = 1.0 / 2.0;
        const double gauss_points_master[gauss_num][3] =
        {
            { 1.0 - 2.0 * gauss_a, gauss_a, gauss_a },
            { 1.0 - gauss_b,       gauss_b, 0.0     },
            { 1.0 - 2.0 * gauss_b, gauss_b, gauss_b },
            { 1.0 - gauss_b,       0.0,     gauss_b },
            { 0.0,                 1.0,     0.0 },
            { 0.0,                 0.0,     1.0 },
            { 1.0,                 0.0,     0.0 }
        };
    }
}

namespace tr_integration_5
{
    namespace tr_integration
    {
        const double gauss_f = 9.0 / 80.0;
        const double gauss_g = (155.0 + sqrt(15.0)) / 2400.0;
        const double gauss_h = (155.0 - sqrt(15.0)) / 2400.0;
        const double gauss_weights[gauss_num] =
        {
            gauss_f,
            gauss_g,
            gauss_g,
            gauss_g,
            gauss_h,
            gauss_h,
            gauss_h
        };
        const double gauss_a = (6.0 + sqrt(15.0)) / 21.0;
        const double gauss_b = (6.0 - sqrt(15.0)) / 21.0;
        const double gauss_c = (9.0 + 2.0 * sqrt(15.0)) / 21.0;
        const double gauss_d = (9.0 - 2.0 * sqrt(15.0)) / 21.0;
        const double gauss_e = 1.0 / 3.0;
        const double gauss_points_master[gauss_num][3] =
        {
            { 1.0 - 2.0 * gauss_e,     gauss_e, gauss_e },
            { 1.0 - gauss_a - gauss_d, gauss_a, gauss_d },
            { 1.0 - 2.0 * gauss_a,     gauss_a, gauss_a },
            { 1.0 - gauss_a - gauss_d, gauss_d, gauss_a },
            { 1.0 - gauss_b - gauss_c, gauss_c, gauss_b },
            { 1.0 - gauss_b - gauss_c, gauss_b, gauss_c },
            { 1.0 - 2.0 * gauss_b,     gauss_b, gauss_b }
        };
    }
}

