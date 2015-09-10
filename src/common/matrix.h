#if !defined MATRIX_H_INCLUDED
#define MATRIX_H_INCLUDED

#include "../common/common.h"

// Класс "статический массив"
template<typename type, size_t dimension = 0>
class array_t
{
public:
    type & operator [] (size_t i)
    {
        return a[i];
    }
    const type & operator [] (size_t i) const
    {
        return a[i];
    }
    array_t()
    {
        for(size_t i = 0; i < dimension; i++)
            a[i] = type();
    }
protected:
    type a[dimension];
};

// Класс "динамический массив"
template<typename type>
class array_t<type, 0>
{
public:
    array_t(size_t dimension)
    {
        this->dimension = dimension;
        a = new type [dimension];
        for(size_t i = 0; i < dimension; i++)
            a[i] = type();
    }
    array_t(const array_t & other)
    {
        dimension = other.dimension;
        a = new type [dimension];
        for(size_t i = 0; i < dimension; i++)
            a[i] = other.a[i];
    }
    ~array_t()
    {
        delete [] a;
    }
    array_t & operator = (const array_t & other)
    {
        if(this->a != other.a)
        {
            delete [] a;
            dimension = other.dimension;
            a = new type [dimension];
            for(size_t i = 0; i < dimension; i++)
                a[i] = other.a[i];
        }
        return * this;
    }
    type & operator [] (size_t i)
    {
        return a[i];
    }
    const type & operator [] (size_t i) const
    {
        return a[i];
    }
protected:
    type * a;
    size_t dimension;
};

// Класс "статическая матрица"
template<typename type, size_t dimension_row = 0, size_t dimension_col = 0>
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
    matrix_t()
    {
        for(size_t i = 0; i < dimension_row; i++)
            for(size_t j = 0; j < dimension_col; j++)
                a[i][j] = type();
    }
protected:
    type a[dimension_row][dimension_col];
};

// Класс "динамическая матрица"
template<typename type>
class matrix_t<type, 0, 0>
{
public:
    matrix_t(size_t dimension_row, size_t dimension_col)
    {
        this->dimension_row = dimension_row;
        this->dimension_col = dimension_col;
        size_t len = dimension_row * dimension_col;
        a = new type [len];
        for(size_t i = 0; i < len; i++)
            a[i] = type();
    }
    matrix_t(const matrix_t & other)
    {
        dimension_row = other.dimension_row;
        dimension_col = other.dimension_col;
        size_t len = dimension_row * dimension_col;
        a = new type [len];
        for(size_t i = 0; i < len; i++)
            a[i] = other.a[i];
    }
    ~matrix_t()
    {
        delete [] a;
    }
    matrix_t & operator = (const matrix_t & other)
    {
        if(this->a != other.a)
        {
            delete [] a;
            dimension_row = other.dimension_row;
            dimension_col = other.dimension_col;
            size_t len = dimension_row * dimension_col;
            a = new type [len];
            for(size_t i = 0; i < len; i++)
                a[i] = other.a[i];
        }
        return * this;
    }
    type * operator [] (size_t i)
    {
        return a + dimension_row * i;
    }
    const type * operator [] (size_t i) const
    {
        return a + dimension_row * i;
    }
protected:
    type * a;
    size_t dimension_row, dimension_col;
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
