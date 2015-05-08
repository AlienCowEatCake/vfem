#ifndef CUBATURES_H_INCLUDED
#define CUBATURES_H_INCLUDED

#include <cstring>
#include <cmath>

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

// Интегрирование особо высоких порядков на тетраэдрах
// http://lsec.cc.ac.cn/~tcui/myinfo/paper/quad.pdf
// http://lsec.cc.ac.cn/phg/download.htm

namespace tet_integration_5
{
    namespace tet_integration
    {
        static const size_t gauss_num = 14;
        extern const double gauss_weights[gauss_num];
        extern const double gauss_points_master[gauss_num][4];
    }
}

namespace tet_integration_5_2 // Stroud T3:5-1 p315
{
    namespace tet_integration
    {
        static const size_t gauss_num = 15;
        extern const double gauss_weights[gauss_num];
        extern const double gauss_points_master[gauss_num][4];
    }
}

namespace tet_integration_6
{
    namespace tet_integration
    {
        static const size_t gauss_num = 24;
        extern const double gauss_weights[gauss_num];
        extern const double gauss_points_master[gauss_num][4];
    }
}

namespace tet_integration_7
{
    namespace tet_integration
    {
        static const size_t gauss_num = 35;
        extern const double gauss_weights[gauss_num];
        extern const double gauss_points_master[gauss_num][4];
    }
}

namespace tet_integration_8
{
    namespace tet_integration
    {
        static const size_t gauss_num = 46;
        extern const double gauss_weights[gauss_num];
        extern const double gauss_points_master[gauss_num][4];
    }
}

namespace tet_integration_9
{
    namespace tet_integration
    {
        static const size_t gauss_num = 59;
        extern const double gauss_weights[gauss_num];
        extern const double gauss_points_master[gauss_num][4];
    }
}

namespace tet_integration_10
{
    namespace tet_integration
    {
        static const size_t gauss_num = 79;
        extern const double gauss_weights[gauss_num];
        extern const double gauss_points_master[gauss_num][4];
    }
}

namespace tet_integration_11
{
    namespace tet_integration
    {
        static const size_t gauss_num = 96;
        extern const double gauss_weights[gauss_num];
        extern const double gauss_points_master[gauss_num][4];
    }
}

namespace tet_integration_12
{
    namespace tet_integration
    {
        static const size_t gauss_num = 127;
        extern const double gauss_weights[gauss_num];
        extern const double gauss_points_master[gauss_num][4];
    }
}

namespace tet_integration_13
{
    namespace tet_integration
    {
        static const size_t gauss_num = 149;
        extern const double gauss_weights[gauss_num];
        extern const double gauss_points_master[gauss_num][4];
    }
}

namespace tet_integration_14
{
    namespace tet_integration
    {
        static const size_t gauss_num = 236;
        extern const double gauss_weights[gauss_num];
        extern const double gauss_points_master[gauss_num][4];
    }
}

// ============================================================================

// Интегрирование особо высоких порядков на треугольниках
// http://lsec.cc.ac.cn/~tcui/myinfo/paper/quad.pdf
// http://lsec.cc.ac.cn/phg/download.htm

namespace tr_integration_6
{
    namespace tr_integration
    {
        static const size_t gauss_num = 12;
        extern const double gauss_weights[gauss_num];
        extern const double gauss_points_master[gauss_num][3];
    }
}

namespace tr_integration_7
{
    namespace tr_integration
    {
        static const size_t gauss_num = 15;
        extern const double gauss_weights[gauss_num];
        extern const double gauss_points_master[gauss_num][3];
    }
}

namespace tr_integration_8
{
    namespace tr_integration
    {
        static const size_t gauss_num = 16;
        extern const double gauss_weights[gauss_num];
        extern const double gauss_points_master[gauss_num][3];
    }
}

namespace tr_integration_9
{
    namespace tr_integration
    {
        static const size_t gauss_num = 19;
        extern const double gauss_weights[gauss_num];
        extern const double gauss_points_master[gauss_num][3];
    }
}

namespace tr_integration_10
{
    namespace tr_integration
    {
        static const size_t gauss_num = 25;
        extern const double gauss_weights[gauss_num];
        extern const double gauss_points_master[gauss_num][3];
    }
}

namespace tr_integration_11
{
    namespace tr_integration
    {
        static const size_t gauss_num = 28;
        extern const double gauss_weights[gauss_num];
        extern const double gauss_points_master[gauss_num][3];
    }
}

namespace tr_integration_12
{
    namespace tr_integration
    {
        static const size_t gauss_num = 33;
        extern const double gauss_weights[gauss_num];
        extern const double gauss_points_master[gauss_num][3];
    }
}

namespace tr_integration_13
{
    namespace tr_integration
    {
        static const size_t gauss_num = 37;
        extern const double gauss_weights[gauss_num];
        extern const double gauss_points_master[gauss_num][3];
    }
}

namespace tr_integration_14
{
    namespace tr_integration
    {
        static const size_t gauss_num = 42;
        extern const double gauss_weights[gauss_num];
        extern const double gauss_points_master[gauss_num][3];
    }
}

// a 49-point rule communicated by Freddie Witherden (freddie@witherden.org)
// 2D 49 point order 15 rule (1:4:6, 346 tries, 3288 evals), error = 4.251e-17
namespace tr_integration_15
{
    namespace tr_integration
    {
        static const size_t gauss_num = 49;
        extern const double gauss_weights[gauss_num];
        extern const double gauss_points_master[gauss_num][3];
    }
}

namespace tr_integration_16
{
    namespace tr_integration
    {
        static const size_t gauss_num = 55;
        extern const double gauss_weights[gauss_num];
        extern const double gauss_points_master[gauss_num][3];
    }
}

// A 60-point order 17 rule found by using the number of points reported in:
// H. Xiao and Z. Gimbutas,
// A numerical algorithm for the construction of efficient quadrature
// rules in two and higher dimensions,
// Computers and Mathematics with Applications, 59 (2010), 663-676
namespace tr_integration_17
{
    namespace tr_integration
    {
        static const size_t gauss_num = 60;
        extern const double gauss_weights[gauss_num];
        extern const double gauss_points_master[gauss_num][3];
    }
}

// A 67-point order 18 rule found by using the number of points reported in:
// H. Xiao and Z. Gimbutas,
// A numerical algorithm for the construction of efficient quadrature
// rules in two and higher dimensions,
// Computers and Mathematics with Applications, 59 (2010), 663-676
namespace tr_integration_18
{
    namespace tr_integration
    {
        static const size_t gauss_num = 67;
        extern const double gauss_weights[gauss_num];
        extern const double gauss_points_master[gauss_num][3];
    }
}

// Note: the rule QUAD_2D_P19 was taken from the book by
// P. Solin, K. Segeth, and I. Dolezel,
// "Higer-order Finite Element Methods",
// Chapman and Hall/CRC Press, 2003.
namespace tr_integration_19
{
    namespace tr_integration
    {
        static const size_t gauss_num = 73;
        extern const double gauss_weights[gauss_num];
        extern const double gauss_points_master[gauss_num][3];
    }
}

// Rules 20-24, 26-29 are found by Simone Weikl (simone.weikl@googlemail.com),
// Diploma Thesis, Technische Universiat Munchen, Zentrum Mathematik, 2011
namespace tr_integration_20
{
    namespace tr_integration
    {
        static const size_t gauss_num = 82;
        extern const double gauss_weights[gauss_num];
        extern const double gauss_points_master[gauss_num][3];
    }
}

namespace tr_integration_21
{
    namespace tr_integration
    {
        static const size_t gauss_num = 90;
        extern const double gauss_weights[gauss_num];
        extern const double gauss_points_master[gauss_num][3];
    }
}

namespace tr_integration_22
{
    namespace tr_integration
    {
        static const size_t gauss_num = 99;
        extern const double gauss_weights[gauss_num];
        extern const double gauss_points_master[gauss_num][3];
    }
}

namespace tr_integration_23
{
    namespace tr_integration
    {
        static const size_t gauss_num = 103;
        extern const double gauss_weights[gauss_num];
        extern const double gauss_points_master[gauss_num][3];
    }
}

namespace tr_integration_24
{
    namespace tr_integration
    {
        static const size_t gauss_num = 117;
        extern const double gauss_weights[gauss_num];
        extern const double gauss_points_master[gauss_num][3];
    }
}

// The rule below is reported in the paper:
// S. Wandzura and H. Xiao, Symmetric quadrature rules on a triangle,
// Computers and Mathematics with Applications, 45 (2003), 1829ЁC1840,
// and was communicated by Don Wilton (dwilton@mindspring.com)
namespace tr_integration_25
{
    namespace tr_integration
    {
        static const size_t gauss_num = 126;
        extern const double gauss_weights[gauss_num];
        extern const double gauss_points_master[gauss_num][3];
    }
}

namespace tr_integration_26
{
    namespace tr_integration
    {
        static const size_t gauss_num = 133;
        extern const double gauss_weights[gauss_num];
        extern const double gauss_points_master[gauss_num][3];
    }
}

namespace tr_integration_27
{
    namespace tr_integration
    {
        static const size_t gauss_num = 145;
        extern const double gauss_weights[gauss_num];
        extern const double gauss_points_master[gauss_num][3];
    }
}

namespace tr_integration_28
{
    namespace tr_integration
    {
        static const size_t gauss_num = 154;
        extern const double gauss_weights[gauss_num];
        extern const double gauss_points_master[gauss_num][3];
    }
}

namespace tr_integration_29
{
    namespace tr_integration
    {
        static const size_t gauss_num = 166;
        extern const double gauss_weights[gauss_num];
        extern const double gauss_points_master[gauss_num][3];
    }
}

// ============================================================================

#endif // CUBATURES_H_INCLUDED
