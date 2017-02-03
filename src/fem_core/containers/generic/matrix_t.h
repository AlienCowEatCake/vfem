#if !defined(CONTAINERS_GENERIC_MATRIX_T_H_INCLUDED)
#define CONTAINERS_GENERIC_MATRIX_T_H_INCLUDED

#include <cstring>
#include <cassert>

namespace fem_core { namespace containers { namespace generic {

// *************************************************************************************************

/**
  * @brief Класс "статическая матрица"
  */
template<typename type, std::size_t dimension_row = 0, std::size_t dimension_col = 0>
class matrix_t
{
    friend class matrix_t<type, 0, 0>;

public:

    /**
      * @brief Конструктор по умолчанию
      */
    matrix_t()
    {
        for(std::size_t i = 0; i < dimension_row; i++)
            for(std::size_t j = 0; j < dimension_col; j++)
                m_values[i][j] = type();
    }

    /**
      * @brief Определение количества строк матрицы
      * @return Количество строк матрицы
      */
    std::size_t size_rows() const
    {
        return dimension_row;
    }

    /**
      * @brief Определение количества столбцов матрицы
      * @return Количество столбцов матрицы
      */
    std::size_t size_cols() const
    {
        return dimension_col;
    }

    /**
      * @brief Оператор [] (константный)
      * @param[in] i Индекс
      * @return Константный указатель на строку матрицы с индексом i
      */
    const type * operator [] (std::size_t i) const
    {
        assert(i < dimension_row);
        return m_values[i];
    }

    /**
      * @brief Оператор []
      * @param[in] i Индекс
      * @return Указатель на строку матрицы с индексом i
      */
    type * operator [] (std::size_t i)
    {
        assert(i < dimension_row);
        return m_values[i];
    }

protected:

    /**
      * @brief Хранимые значения, двумерный статический массив на стеке
      */
    type m_values[dimension_row][dimension_col];
};

// *************************************************************************************************

/**
  * @brief Класс "динамическая матрица"
  */
template<typename type>
class matrix_t<type, 0, 0>
{
public:

    /**
      * @brief Конструктор с опциональным размером
      */
    matrix_t(std::size_t dimension_row = 0, std::size_t dimension_col = 0)
        : m_values(NULL)
    {
        resize(dimension_row, dimension_col);
    }

    /**
      * @brief Конструктор копирования
      * @param[in] other Копируемая матрица
      */
    matrix_t(const matrix_t & other)
        : m_values(NULL)
    {
        copy_from(other);
    }

    /**
      * @brief Конструктор копирования из неявно приводимой матрицы
      * @param[in] other Копируемая неявно приводимая матрица
      */
    template<typename type_other, std::size_t dimension_row_other, std::size_t dimension_col_other>
    matrix_t(const matrix_t<type_other, dimension_row_other, dimension_col_other> & other)
        : m_values(NULL)
    {
        copy_from(other);
    }

    /**
      * @brief Деструктор
      */
    ~matrix_t()
    {
        delete [] m_values;
    }

    /**
      * @brief Оператор присваивания
      * @param[in] other Матрица, которую нужно присвоить
      * @return *this
      */
    matrix_t & operator = (const matrix_t & other)
    {
        if(m_values != other.m_values)
            copy_from(other);
        return * this;
    }

    /**
      * @brief Оператор присваивания неявно приводимой динамической матрицы
      * @param[in] other Неявно приводимая динамическая матрица, которую нужно присвоить
      * @return *this
      */
    template<typename type_other>
    matrix_t & operator = (const matrix_t<type_other, 0, 0> & other)
    {
        if(m_values != other.m_values)
            copy_from(other);
        return * this;
    }

    /**
      * @brief Оператор присваивания неявно приводимой статической матрицы
      * @param[in] other Неявно приводимая статическая матрица, которую нужно присвоить
      * @return *this
      */
    template<typename type_other, std::size_t dimension_row_other, std::size_t dimension_col_other>
    matrix_t & operator = (const matrix_t<type_other, dimension_row_other, dimension_col_other> & other)
    {
        copy_from(other);
        return * this;
    }

    /**
      * @brief Определение количества строк матрицы
      * @return Количество строк матрицы
      */
    std::size_t size_rows() const
    {
        return m_dimension_row;
    }

    /**
      * @brief Определение количества столбцов матрицы
      * @return Количество столбцов матрицы
      */
    std::size_t size_cols() const
    {
        return m_dimension_col;
    }

    /**
      * @brief Изменение размера матрицы (с потерей данных)
      * @param dimension_row Новое количество строк
      * @param dimension_col Новое количество столбцов
      */
    void resize(std::size_t dimension_row, std::size_t dimension_col)
    {
        realloc(dimension_row, dimension_col);
        const std::size_t len = dimension_row * dimension_col;
        for(std::size_t i = 0; i < len; i++)
            m_values[i] = type();
    }

    /**
      * @brief Оператор [] (константный)
      * @param[in] i Индекс
      * @return Константный указатель на строку матрицы с индексом i
      */
    const type * operator [] (std::size_t i) const
    {
        assert(i < m_dimension_row);
        return m_values + m_dimension_col * i;
    }

    /**
      * @brief Оператор []
      * @param[in] i Индекс
      * @return Указатель на строку матрицы с индексом i
      */
    type * operator [] (std::size_t i)
    {
        assert(i < m_dimension_row);
        return m_values + m_dimension_col * i;
    }

protected:

    /**
      * @brief Хранимые значения, динамический массив
      */
    type * m_values;
    /**
      * @brief Количество строк матрицы
      */
    std::size_t m_dimension_row;
    /**
      * @brief Количество столбцов матрицы
      */
    std::size_t m_dimension_col;

    /**
      * @brief Изменение размера матрицы (с потерей данных)
      * @param dimension_row Новое количество строк
      * @param dimension_col Новое количество столбцов
      */
    void realloc(std::size_t dimension_row, std::size_t dimension_col)
    {
        m_dimension_row = dimension_row;
        m_dimension_col = dimension_col;
        std::size_t len = dimension_row * dimension_col;
        delete [] m_values;
        if(len == 0)
            m_values = NULL;
        else
            m_values = new type [len];
    }

    /**
      * @brief Копирование из другой неявно приводимой динамической матрицы
      * @param[in] other другая неявно приводимая динамическая матрица
      */
    template<typename type_other>
    void copy_from(const matrix_t<type_other, 0, 0> & other)
    {
        realloc(other.m_dimension_row, other.m_dimension_col);
        const std::size_t len = m_dimension_row * m_dimension_col;
        for(std::size_t i = 0; i < len; i++)
            m_values[i] = other.m_values[i];
    }

    /**
      * @brief Копирование из другой неявно приводимой статической матрицы
      * @param[in] other другая неявно приводимая статическая матрица
      */
    template<typename type_other, std::size_t dimension_row_other, std::size_t dimension_col_other>
    void copy_from(const matrix_t<type_other, dimension_row_other, dimension_col_other> & other)
    {
        realloc(dimension_row_other, dimension_col_other);
        for(std::size_t i = 0; i < m_dimension_row; i++)
        {
            type * str = m_values + m_dimension_col * i;
            for(std::size_t j = 0; j < m_dimension_col; j++)
                str[j] = other.m_values[i][j];
        }
    }
};

// *************************************************************************************************

/**
  * @brief Определитель матрицы 2x2
  * @param[in] A Исходная матрица
  * @return Определитель исходной матрицы
  */
template<typename type>
type determenant(const matrix_t<type, 2, 2> & A)
{
    return A[0][0] * A[1][1] -  A[0][1] * A[1][0];
}

/**
  * @brief Определитель матрицы 3x3
  * @param[in] A Исходная матрица
  * @return Определитель исходной матрицы
  */
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

/**
  * @brief Определитель матрицы 4x4
  * @param[in] A Исходная матрица
  * @return Определитель исходной матрицы
  */
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

// *************************************************************************************************

/**
  * @brief Обращение матрицы 2x2
  * @param[in] A Исходная матрица
  * @param[out] detA Определитель исходной матрицы
  * @return Матрица, обратная к исходной
  */
template<typename type>
matrix_t<type, 2, 2> inverse(const matrix_t<type, 2, 2> & A, type & detA)
{
    matrix_t<type, 2, 2> Inv_A;
    detA = determenant(A);

    Inv_A[0][0] =   A[1][1] / detA;
    Inv_A[0][1] = - A[0][1] / detA;
    Inv_A[1][0] = - A[1][0] / detA;
    Inv_A[1][1] =   A[0][0] / detA;

    return Inv_A;
}

/**
  * @brief Обращение матрицы 3x3
  * @param[in] A Исходная матрица
  * @param[out] detA Определитель исходной матрицы
  * @return Матрица, обратная к исходной
  */
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

/**
  * @brief Обращение матрицы 4x4
  * @param[in] A Исходная матрица
  * @param[out] detA Определитель исходной матрицы
  * @return Матрица, обратная к исходной
  */
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

// *************************************************************************************************

}}} // namespace fem_core::containers::generic

#endif // CONTAINERS_GENERIC_MATRIX_T_H_INCLUDED
