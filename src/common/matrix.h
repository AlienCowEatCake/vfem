#if !defined MATRIX_H_INCLUDED
#define MATRIX_H_INCLUDED

#include "../common/common.h"

// Класс "статический массив"
template<typename type, size_t dimension = 0>
class array_t
{
    friend class array_t<type, 0>;

protected:

    type a[dimension];

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
};


// Класс "динамический массив"
template<typename type>
class array_t<type, 0>
{
protected:

    type * a;
    size_t dimension;

    void realloc(size_t dimension)
    {
        this->dimension = dimension;
        delete [] a;
        if(dimension == 0)
            a = NULL;
        else
            a = new type [dimension];
    }

    template<typename type_other>
    void copy_from(const array_t<type_other, 0> & other)
    {
        realloc(other.dimension);
        for(size_t i = 0; i < dimension; i++)
            a[i] = other.a[i];
    }

    template<typename type_other, size_t dimension_other>
    void copy_from(const array_t<type_other, dimension_other> & other)
    {
        realloc(dimension_other);
        for(size_t i = 0; i < dimension; i++)
            a[i] = other.a[i];
    }

public:

    void resize(size_t dimension)
    {
        realloc(dimension);
        for(size_t i = 0; i < dimension; i++)
            a[i] = type();
    }

    array_t(size_t dimension = 0)
    {
        a = NULL;
        resize(dimension);
    }

    array_t(const array_t & other)
    {
        a = NULL;
        copy_from(other);
    }

    template<typename type_other, size_t dimension_other>
    array_t(const array_t<type_other, dimension_other> & other)
    {
        a = NULL;
        copy_from(other);
    }

    ~array_t()
    {
        delete [] a;
    }

    array_t & operator = (const array_t & other)
    {
        if(this->a != other.a)
            copy_from(other);
        return * this;
    }

    template<typename type_other>
    array_t & operator = (const array_t<type_other, 0> & other)
    {
        if(this->a != other.a)
            copy_from(other);
        return * this;
    }

    template<typename type_other, size_t dimension_other>
    array_t & operator = (const array_t<type_other, dimension_other> & other)
    {
        copy_from(other);
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
};


// Класс "статическая матрица"
template<typename type, size_t dimension_row = 0, size_t dimension_col = 0>
class matrix_t
{
    friend class matrix_t<type, 0, 0>;

protected:

    type a[dimension_row][dimension_col];

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
};


// Класс "динамическая матрица"
template<typename type>
class matrix_t<type, 0, 0>
{
protected:

    type * a;
    size_t dimension_row, dimension_col;

    void realloc(size_t dimension_row, size_t dimension_col)
    {
        this->dimension_row = dimension_row;
        this->dimension_col = dimension_col;
        size_t len = dimension_row * dimension_col;
        delete [] a;
        if(len == 0)
            a = NULL;
        else
            a = new type [len];
    }

    template<typename type_other>
    void copy_from(const matrix_t<type_other, 0, 0> & other)
    {
        realloc(other.dimension_row, other.dimension_col);
        size_t len = dimension_row * dimension_col;
        for(size_t i = 0; i < len; i++)
            a[i] = other.a[i];
    }

    template<typename type_other, size_t dimension_row_other, size_t dimension_col_other>
    void copy_from(const matrix_t<type_other, dimension_row_other, dimension_col_other> & other)
    {
        realloc(dimension_row_other, dimension_col_other);
        for(size_t i = 0; i < dimension_row; i++)
        {
            type * str = a + dimension_col * i;
            for(size_t j = 0; j < dimension_col; j++)
                str[j] = other.a[i][j];
        }
    }

public:

    void resize(size_t dimension_row, size_t dimension_col)
    {
        realloc(dimension_row, dimension_col);
        size_t len = dimension_row * dimension_col;
        for(size_t i = 0; i < len; i++)
            a[i] = type();
    }

    matrix_t(size_t dimension_row = 0, size_t dimension_col = 0)
    {
        a = NULL;
        resize(dimension_row, dimension_col);
    }

    matrix_t(const matrix_t & other)
    {
        a = NULL;
        copy_from(other);
    }

    template<typename type_other, size_t dimension_row_other, size_t dimension_col_other>
    matrix_t(const matrix_t<type_other, dimension_row_other, dimension_col_other> & other)
    {
        a = NULL;
        copy_from(other);
    }

    ~matrix_t()
    {
        delete [] a;
    }

    matrix_t & operator = (const matrix_t & other)
    {
        if(this->a != other.a)
            copy_from(other);
        return * this;
    }

    template<typename type_other>
    matrix_t & operator = (const matrix_t<type_other, 0, 0> & other)
    {
        if(this->a != other.a)
            copy_from(other);
        return * this;
    }

    template<typename type_other, size_t dimension_row_other, size_t dimension_col_other>
    matrix_t & operator = (const matrix_t<type_other, dimension_row_other, dimension_col_other> & other)
    {
        copy_from(other);
        return * this;
    }

    type * operator [] (size_t i)
    {
        return a + dimension_col * i;
    }

    const type * operator [] (size_t i) const
    {
        return a + dimension_col * i;
    }
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
