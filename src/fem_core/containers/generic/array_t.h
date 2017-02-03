#if !defined(CONTAINERS_GENERIC_ARRAY_T_H_INCLUDED)
#define CONTAINERS_GENERIC_ARRAY_T_H_INCLUDED

#include <cstring>
#include <cassert>

namespace fem_core { namespace containers { namespace generic {

// *************************************************************************************************

/**
  * @brief Класс "статический массив"
  */
template<typename type, std::size_t dimension = 0>
class array_t
{
    friend class array_t<type, 0>;

public:

    /**
      * @brief Конструктор по умолчанию
      */
    array_t()
    {
        for(std::size_t i = 0; i < dimension; i++)
            m_values[i] = type();
    }

    /**
      * @brief Определение размера массива
      * @return Размер массива
      */
    std::size_t size() const
    {
        return dimension;
    }

    /**
      * @brief Оператор [] (константный)
      * @param[in] i Индекс
      * @return Константная ссылка на элемент массива с индексом i
      */
    const type & operator [] (std::size_t i) const
    {
        assert(i < dimension);
        return m_values[i];
    }

    /**
      * @brief Оператор []
      * @param[in] i Индекс
      * @return Ссылка на элемент массива с индексом i
      */
    type & operator [] (std::size_t i)
    {
        assert(i < dimension);
        return m_values[i];
    }

    /**
      * @brief Получение констаниного сырого указателя на данные
      * @return Константный сырой указатель на данные
      */
    const type * data() const
    {
        return m_values;
    }

    /**
      * @brief Получение сырого указателя на данные
      * @return Сырой указатель на данные
      */
    type * data()
    {
        return m_values;
    }

protected:

    /**
      * @brief Хранимые значения, статический массив на стеке
      */
    type m_values[dimension];
};

// *************************************************************************************************

/**
  * @brief Класс "динамический массив"
  */
template<typename type>
class array_t<type, 0>
{
public:

    /**
      * @brief Конструктор с опциональным размером
      */
    array_t(std::size_t dimension = 0)
        : m_values(NULL)
    {
        resize(dimension);
    }

    /**
      * @brief Конструктор копирования
      * @param[in] other Копируемый массив
      */
    array_t(const array_t & other)
        : m_values(NULL)
    {
        copy_from(other);
    }

    /**
      * @brief Конструктор копирования из неявно приводимого массива
      * @param[in] other Копируемый неявно приводимый массив
      */
    template<typename type_other, std::size_t dimension_other>
    array_t(const array_t<type_other, dimension_other> & other)
        : m_values(NULL)
    {
        copy_from(other);
    }

    /**
      * @brief Деструктор
      */
    ~array_t()
    {
        delete [] m_values;
    }

    /**
      * @brief Оператор присваивания
      * @param[in] other Массив, который нужно присвоить
      * @return *this
      */
    array_t & operator = (const array_t & other)
    {
        if(m_values != other.m_values)
            copy_from(other);
        return * this;
    }

    /**
      * @brief Оператор присваивания неявно приводимого динамического массива
      * @param[in] other Неявно приводимый динамический массив, который нужно присвоить
      * @return *this
      */
    template<typename type_other>
    array_t & operator = (const array_t<type_other, 0> & other)
    {
        if(m_values != other.m_values)
            copy_from(other);
        return * this;
    }

    /**
      * @brief Оператор присваивания неявно приводимого статического массива
      * @param[in] other Неявно приводимый статический массив, который нужно присвоить
      * @return *this
      */
    template<typename type_other, std::size_t dimension_other>
    array_t & operator = (const array_t<type_other, dimension_other> & other)
    {
        copy_from(other);
        return * this;
    }

    /**
      * @brief Определение размера массива
      * @return Размер массива
      */
    std::size_t size()
    {
        return m_dimension;
    }

    /**
      * @brief Изменение размера массива (с потерей данных)
      * @param[in] dimension Новый размер массива
      */
    void resize(std::size_t dimension)
    {
        realloc(dimension);
        for(std::size_t i = 0; i < dimension; i++)
            m_values[i] = type();
    }

    /**
      * @brief Оператор [] (константный)
      * @param[in] i Индекс
      * @return Константная ссылка на элемент массива с индексом i
      */
    const type & operator [] (std::size_t i) const
    {
        assert(i < m_dimension);
        return m_values[i];
    }

    /**
      * @brief Оператор []
      * @param[in] i Индекс
      * @return Ссылка на элемент массива с индексом i
      */
    type & operator [] (std::size_t i)
    {
        assert(i < m_dimension);
        return m_values[i];
    }

    /**
      * @brief Получение констаниного сырого указателя на данные
      * @return Константный сырой указатель на данные
      */
    const type * data() const
    {
        return m_values;
    }

    /**
      * @brief Получение сырого указателя на данные
      * @return Сырой указатель на данные
      */
    type * data()
    {
        return m_values;
    }

protected:

    /**
      * @brief Хранимые значения, динамический массив
      */
    type * m_values;
    /**
      * @brief Размерность массива хранимых значений
      */
    std::size_t m_dimension;

    /**
      * @brief Копирование из другого неявно приводимого динамического массива
      * @param[in] other другой неявно приводимый динамический массив
      */
    template<typename type_other>
    void copy_from(const array_t<type_other, 0> & other)
    {
        realloc(other.m_dimension);
        for(std::size_t i = 0; i < m_dimension; i++)
            m_values[i] = other.m_values[i];
    }

    /**
      * @brief Копирование из другого неявно приводимого статического массива
      * @param[in] other другой неявно приводимый статический массив
      */
    template<typename type_other, std::size_t dimension_other>
    void copy_from(const array_t<type_other, dimension_other> & other)
    {
        realloc(dimension_other);
        for(std::size_t i = 0; i < m_dimension; i++)
            m_values[i] = other.m_values[i];
    }

    /**
      * @brief Изменение размера массива (с потерей данных)
      * @param[in] dimension Новый размер массива
      */
    void realloc(std::size_t dimension)
    {
        m_dimension = dimension;
        delete [] m_values;
        if(dimension == 0)
            m_values = NULL;
        else
            m_values = new type [dimension];
    }
};

// *************************************************************************************************

}}} // namespace fem_core::containers::generic

#endif // CONTAINERS_GENERIC_ARRAY_T_H_INCLUDED
