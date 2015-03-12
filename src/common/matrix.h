#if !defined MATRIX_H_INCLUDED
#define MATRIX_H_INCLUDED

#include "../common/common.h"

// Немного typedef'ов
template<typename type, size_t dimension> class array_t;
template<typename type, size_t dimension_row, size_t dimension_col> class matrix_t;
typedef matrix_t<double, 3, 3> matrix3;
typedef matrix_t<double, 4, 4> matrix4;
// # DELETEME
typedef matrix_t<double,  6,  6> matrix6;
typedef matrix_t<double, 12, 12> matrix12;
typedef array_t<double,  6> array6;
typedef array_t<double, 12> array12;
// # EODELETEME
typedef matrix_t<complex<double>, 3, 3> cmatrix3;
typedef matrix_t<complex<double>, 4, 4> cmatrix4;
// # DELETEME
typedef matrix_t<complex<double>,  6, 6> cmatrix6;
typedef matrix_t<complex<double>, 12, 12> cmatrix12;
typedef array_t<complex<double>,  6> carray6;
typedef array_t<complex<double>, 12> carray12;
// # EODELETEME

// Класс "статический массив"
template<typename type, size_t dimension>
class array_t
{
public:
    type & operator [] (size_t i)
    {
        return a[i];
    }
    type operator [] (size_t i) const
    {
        return a[i];
    }
protected:
    type a[dimension];
};

// Класс "статическая матрица"
template<typename type, size_t dimension_row, size_t dimension_col>
class matrix_t
{
public:
    type * operator [] (size_t i)
    {
        return a[i];
    }
    const type * operator [] (size_t i) const
    {
        return a[i];
    }
protected:
    type a[dimension_row][dimension_col];
};

// Определитель матрицы 3x3
template<typename type>
type determenant(const matrix_t<type, 3, 3> & A)
{
    return A[0][0] * A[1][1] * A[2][2] +
           A[1][0] * A[2][1] * A[0][2] +
           A[0][1] * A[1][2] * A[2][0] -
           A[2][0] * A[1][1] * A[0][2] -
           A[2][1] * A[1][2] * A[0][0] -
           A[1][0] * A[0][1] * A[2][2];
}

// Определитель матрицы 4x4
template<typename type>
type determenant(const matrix_t<type, 4, 4> & A)
{
    return A[2][3] * (A[3][1] * (A[0][0] * A[1][2] - A[0][2] * A[1][0]) +
                      A[0][1] * (A[1][0] * A[3][2] - A[1][2] * A[3][0])) +
           A[2][1] * (A[3][2] * (A[0][0] * A[1][3] - A[0][3] * A[1][0]) +
                      A[3][0] * (A[0][3] * A[1][2] - A[0][2] * A[1][3]) +
                      A[3][3] * (A[0][2] * A[1][0] - A[0][0] * A[1][2])) +
           A[1][1] * (A[2][3] * (A[0][2] * A[3][0] - A[0][0] * A[3][2]) +
                      A[2][2] * (A[0][0] * A[3][3] - A[0][3] * A[3][0])) +
           A[2][2] * A[3][1] * (A[0][3] * A[1][0] - A[0][0] * A[1][3]) +
           A[2][0] * (A[0][1] * (A[1][2] * A[3][3] - A[1][3] * A[3][2]) +
                      A[1][1] * (A[0][3] * A[3][2] - A[0][2] * A[3][3]) +
                      A[3][1] * (A[0][2] * A[1][3] - A[0][3] * A[1][2])) +
           A[0][1] * A[2][2] * (A[1][3] * A[3][0] - A[1][0] * A[3][3]);
}

// Обращение матрицы 3x3
template<typename type>
matrix_t<type, 3, 3> inverse(const matrix_t<type, 3, 3> & A, type & detA)
{
    matrix_t<type, 3, 3> Inv_A;
    detA = determenant(A);

    Inv_A[0][0] = (A[1][1] * A[2][2] - A[2][1] * A[1][2]) / detA;
    Inv_A[0][1] = - (A[0][1] * A[2][2] - A[2][1] * A[0][2]) / detA;
    Inv_A[0][2] = (A[0][1] * A[1][2] - A[1][1] * A[0][2]) / detA;

    Inv_A[1][0] = - (A[1][0] * A[2][2] - A[2][0] * A[1][2]) / detA;
    Inv_A[1][1] = (A[0][0] * A[2][2] - A[2][0] * A[0][2]) / detA;
    Inv_A[1][2] = - (A[0][0] * A[1][2] - A[1][0] * A[0][2]) / detA;

    Inv_A[2][0] = (A[1][0] * A[2][1] - A[2][0] * A[1][1]) / detA;
    Inv_A[2][1] = - (A[0][0] * A[2][1] - A[2][0] * A[0][1]) / detA;
    Inv_A[2][2] = (A[0][0] * A[1][1] - A[1][0] * A[0][1]) / detA;

    return Inv_A;
}

// Обращение матрицы 4x4
template<typename type>
matrix_t<type, 4, 4> inverse(const matrix_t<type, 4, 4> & A, type & detA)
{
    matrix_t<type, 4, 4> Inv_A;
    detA = determenant(A);

    Inv_A[0][0] = (A[3][1] * (A[1][2] * A[2][3] - A[1][3] * A[2][2]) +
                   A[3][2] * (A[1][3] * A[2][1] - A[1][1] * A[2][3]) +
                   A[3][3] * (A[1][1] * A[2][2] - A[1][2] * A[2][1])) / detA;
    Inv_A[0][1] = (A[3][1] * (A[0][3] * A[2][2] - A[0][2] * A[2][3]) +
                   A[3][2] * (A[0][1] * A[2][3] - A[0][3] * A[2][1]) +
                   A[3][3] * (A[0][2] * A[2][1] - A[0][1] * A[2][2])) / detA;
    Inv_A[0][2] = (A[3][1] * (A[0][2] * A[1][3] - A[0][3] * A[1][2]) +
                   A[3][2] * (A[0][3] * A[1][1] - A[0][1] * A[1][3]) +
                   A[3][3] * (A[0][1] * A[1][2] - A[0][2] * A[1][1])) / detA;
    Inv_A[0][3] = (A[2][1] * (A[0][3] * A[1][2] - A[0][2] * A[1][3]) +
                   A[2][2] * (A[0][1] * A[1][3] - A[0][3] * A[1][1]) +
                   A[2][3] * (A[0][2] * A[1][1] - A[0][1] * A[1][2])) / detA;

    Inv_A[1][0] = (A[3][0] * (A[1][3] * A[2][2] - A[1][2] * A[2][3]) +
                   A[3][2] * (A[1][0] * A[2][3] - A[1][3] * A[2][0]) +
                   A[3][3] * (A[1][2] * A[2][0] - A[1][0] * A[2][2])) / detA;
    Inv_A[1][1] = (A[3][0] * (A[0][2] * A[2][3] - A[0][3] * A[2][2]) +
                   A[3][2] * (A[0][3] * A[2][0] - A[0][0] * A[2][3]) +
                   A[3][3] * (A[0][0] * A[2][2] - A[0][2] * A[2][0])) / detA;
    Inv_A[1][2] = (A[3][0] * (A[0][3] * A[1][2] - A[0][2] * A[1][3]) +
                   A[3][2] * (A[0][0] * A[1][3] - A[0][3] * A[1][0]) +
                   A[3][3] * (A[0][2] * A[1][0] - A[0][0] * A[1][2])) / detA;
    Inv_A[1][3] = (A[2][0] * (A[0][2] * A[1][3] - A[0][3] * A[1][2]) +
                   A[2][2] * (A[0][3] * A[1][0] - A[0][0] * A[1][3]) +
                   A[2][3] * (A[0][0] * A[1][2] - A[0][2] * A[1][0])) / detA;

    Inv_A[2][0] = (A[3][0] * (A[1][1] * A[2][3] - A[1][3] * A[2][1]) +
                   A[3][1] * (A[1][3] * A[2][0] - A[1][0] * A[2][3]) +
                   A[3][3] * (A[1][0] * A[2][1] - A[1][1] * A[2][0])) / detA;
    Inv_A[2][1] = (A[3][0] * (A[0][3] * A[2][1] - A[0][1] * A[2][3]) +
                   A[3][1] * (A[0][0] * A[2][3] - A[0][3] * A[2][0]) +
                   A[3][3] * (A[0][1] * A[2][0] - A[0][0] * A[2][1])) / detA;
    Inv_A[2][2] = (A[3][0] * (A[0][1] * A[1][3] - A[0][3] * A[1][1]) +
                   A[3][1] * (A[0][3] * A[1][0] - A[0][0] * A[1][3]) +
                   A[3][3] * (A[0][0] * A[1][1] - A[0][1] * A[1][0])) / detA;
    Inv_A[2][3] = (A[2][0] * (A[0][3] * A[1][1] - A[0][1] * A[1][3]) +
                   A[2][1] * (A[0][0] * A[1][3] - A[0][3] * A[1][0]) +
                   A[2][3] * (A[0][1] * A[1][0] - A[0][0] * A[1][1])) / detA;

    Inv_A[3][0] = (A[3][0] * (A[1][2] * A[2][1] - A[1][1] * A[2][2]) +
                   A[3][1] * (A[1][0] * A[2][2] - A[1][2] * A[2][0]) +
                   A[3][2] * (A[1][1] * A[2][0] - A[1][0] * A[2][1])) / detA;
    Inv_A[3][1] = (A[3][0] * (A[0][1] * A[2][2] - A[0][2] * A[2][1]) +
                   A[3][1] * (A[0][2] * A[2][0] - A[0][0] * A[2][2]) +
                   A[3][2] * (A[0][0] * A[2][1] - A[0][1] * A[2][0])) / detA;
    Inv_A[3][2] = (A[3][0] * (A[0][2] * A[1][1] - A[0][1] * A[1][2]) +
                   A[3][1] * (A[0][0] * A[1][2] - A[0][2] * A[1][0]) +
                   A[3][2] * (A[0][1] * A[1][0] - A[0][0] * A[1][1])) / detA;
    Inv_A[3][3] = (A[2][0] * (A[0][1] * A[1][2] - A[0][2] * A[1][1]) +
                   A[2][1] * (A[0][2] * A[1][0] - A[0][0] * A[1][2]) +
                   A[2][2] * (A[0][0] * A[1][1] - A[0][1] * A[1][0])) / detA;

    return Inv_A;
}

#endif // MATRIX_H_INCLUDED
