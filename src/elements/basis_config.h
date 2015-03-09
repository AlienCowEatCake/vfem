#ifndef BASIS_H_INCLUDED
#define BASIS_H_INCLUDED

#include <cstring>
#include <cmath>

#define BASIS_ORDER 1
#define BASIS_FULL  1

// ====== Базис первого неполного порядка =====================================

#if BASIS_ORDER == 1 && BASIS_FULL == 0
namespace basis
{
    static const size_t tet_bf_num = 6;
    static const size_t tr_bf_num = 3;
}
#endif

// ====== Базис первого полного порядка =======================================

#if BASIS_ORDER == 1 && BASIS_FULL == 1
namespace basis
{
    static const size_t tet_bf_num = 12;
    static const size_t tr_bf_num = 6;
}
#endif

// ====== Базис второго неполного порядка =====================================

#if BASIS_ORDER == 2 && BASIS_FULL == 0

#endif

// ====== Базис второго полного порядка =======================================

#if BASIS_ORDER == 2 && BASIS_FULL == 1

#endif

// ============================================================================

#if BASIS_ORDER == 1
#define tet_integration_2 tet_integration
#define tr_integration_2 tr_integration
#endif

#if BASIS_ORDER == 2
#define tet_integration_4 tet_integration
#define tr_integration_5 tr_integration
#endif

// ============================================================================

// Численное интегрирование на тетраэдрах
// https://people.fh-landshut.de/~maurer/numeth/node148.html

namespace tet_integration_2
{
    static const size_t gauss_num = 4;
    static const double gauss_weights[gauss_num] =
    {
        1.0 / 24.0,
        1.0 / 24.0,
        1.0 / 24.0,
        1.0 / 24.0
    };
    static const double gauss_a = (5.0 - sqrt(5.0)) / 20.0;
    static const double gauss_b = (5.0 + 3.0 * sqrt(5.0)) / 20.0;
    static const double gauss_points[gauss_num][4] =
    {
        { 1.0 - gauss_b - 2.0 * gauss_a, gauss_b, gauss_a, gauss_a },
        { 1.0 - gauss_b - 2.0 * gauss_a, gauss_a, gauss_b, gauss_a },
        { 1.0 - gauss_b - 2.0 * gauss_a, gauss_a, gauss_a, gauss_b },
        { 1.0 - 3.0 * gauss_a,           gauss_a, gauss_a, gauss_a }
    };
}

namespace tet_integration_3
{
    static const size_t gauss_num = 5;
    static const double gauss_weights[gauss_num] =
    {
        - 2.0 / 15.0,
        3.0 / 40.0,
        3.0 / 40.0,
        3.0 / 40.0,
        3.0 / 40.0
    };
    static const double gauss_a = 1.0 / 4.0;
    static const double gauss_b = 1.0 / 2.0;
    static const double gauss_c = 1.0 / 6.0;
    static const double gauss_points[gauss_num][4] =
    {
        { 1.0 - 3.0 * gauss_a,           gauss_a, gauss_a, gauss_a },
        { 1.0 - gauss_b - 2.0 * gauss_c, gauss_b, gauss_c, gauss_c },
        { 1.0 - gauss_b - 2.0 * gauss_c, gauss_c, gauss_b, gauss_c },
        { 1.0 - gauss_b - 2.0 * gauss_c, gauss_c, gauss_c, gauss_b },
        { 1.0 - 3.0 * gauss_c,           gauss_c, gauss_c, gauss_c }
    };
}

namespace tet_integration_4
{
    static const size_t gauss_num = 11;
    static const double gauss_weights[gauss_num] =
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
    static const double gauss_a = 1.0 / 4.0;
    static const double gauss_b = 11.0 / 14.0;
    static const double gauss_c = 5.0 / 70.0;
    static const double gauss_d = (1.0 + sqrt(5.0 / 14.0)) / 4.0;
    static const double gauss_e = (1.0 - sqrt(5.0 / 14.0)) / 4.0;
    static const double gauss_points[gauss_num][4] =
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

// ============================================================================

// Численное интегрирование на треугольниках
// https://people.fh-landshut.de/~maurer/numeth/node147.html

namespace tr_integration_2
{
    static const size_t gauss_num = 3;
    static const double gauss_weights[gauss_num] =
    {
        1.0 / 6.0,
        1.0 / 6.0,
        1.0 / 6.0
    };
    static const double gauss_points[gauss_num][3] =
    {
        { 1.0 / 2.0, 1.0 / 2.0, 0.0       },
        { 0.0,       1.0 / 2.0, 1.0 / 2.0 },
        { 1.0 / 2.0, 0.0,       1.0 / 2.0 }
    };
}

namespace tr_integration_3
{
    static const size_t gauss_num = 4;
    static const double gauss_weights[gauss_num] =
    {
        - 9.0 / 32.0,
        25.0 / 96.0,
        25.0 / 96.0,
        25.0 / 96.0
    };
    static const double gauss_a = 1.0 / 3.0;
    static const double gauss_b = 3.0 / 5.0;
    static const double gauss_c = 1.0 / 5.0;
    static const double gauss_points[gauss_num][3] =
    {
        { 1.0 - 2.0 * gauss_a,     gauss_a, gauss_a },
        { 1.0 - gauss_b - gauss_c, gauss_b, gauss_c },
        { 1.0 - gauss_b - gauss_c, gauss_c, gauss_b },
        { 1.0 - 2.0 * gauss_c,     gauss_c, gauss_c }
    };
}

namespace tr_integration_5
{
    static const size_t gauss_num = 7;
    static const double gauss_f = 9.0 / 80.0;
    static const double gauss_g = (155.0 + sqrt(15.0)) / 2400.0;
    static const double gauss_h = (155.0 - sqrt(15.0)) / 2400.0;
    static const double gauss_weights[gauss_num] =
    {
        gauss_f,
        gauss_g,
        gauss_g,
        gauss_g,
        gauss_h,
        gauss_h,
        gauss_h
    };
    static const double gauss_a = (6.0 + sqrt(15.0)) / 21.0;
    static const double gauss_b = (6.0 - sqrt(15.0)) / 21.0;
    static const double gauss_c = (9.0 + 2.0 * sqrt(15.0)) / 21.0;
    static const double gauss_d = (9.0 - 2.0 * sqrt(15.0)) / 21.0;
    static const double gauss_e = 1.0 / 3.0;
    static const double gauss_points[gauss_num][3] =
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

// ============================================================================

#if defined tet_integration_2
#undef tet_integration_2
#endif

#if defined tet_integration_3
#undef tet_integration_3
#endif

#if defined tet_integration_4
#undef tet_integration_4
#endif

#if defined tr_integration_2
#undef tr_integration_2
#endif

#if defined tr_integration_3
#undef tr_integration_3
#endif

#if defined tr_integration_5
#undef tr_integration_5
#endif

#endif // BASIS_H_INCLUDED
