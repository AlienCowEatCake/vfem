#if !defined(CONTAINERS_GEOMETRY_POINT3_T_H_INCLUDED)
#define CONTAINERS_GEOMETRY_POINT3_T_H_INCLUDED

#include <cstring>
#include <cassert>
#include <fstream>

namespace core { namespace containers { namespace geometry {

// *************************************************************************************************

// Структура трехмерная точка
template<typename T>
struct point3_t
{
    // Координаты
    T x, y, z;
    // Номер точки
    std::size_t num;
    // Конструктор по-умолчанию
    point3_t();
    // Конструктор по трем координатам и номеру
    template<typename U>
    point3_t(U n_x, U n_y, U n_z, std::size_t n_num = 0);
    // Конструктор из другой точки
    template<typename U>
    point3_t(const point3_t<U> & p);
    // Операторы типа "скобка"
    T & operator [] (std::size_t i);
    const T & operator [] (std::size_t i) const;
    // Оператор меньше (по номеру)
    template<typename U>
    bool operator < (const point3_t<U> & t) const;
    // Оператор равенства (по номеру)
    template<typename U>
    bool operator == (const point3_t<U> & t) const;
    // Оператор присваивания
    point3_t<T> & operator = (const point3_t<T> & other);
    template<typename U>
    point3_t<T> & operator = (const point3_t<U> & other);
    // Вывод
    template<typename U>
    friend std::ostream & operator << (std::ostream & os, const point3_t<U> & a);
    // Проверка, лежит ли точка внутри параллелепипеда (для дерева)
    bool inside(double x0, double x1, double y0, double y1, double z0, double z1) const;
};

// *************************************************************************************************

// Конструктор по-умолчанию
template<typename T>
point3_t<T>::point3_t()
    : x(static_cast<T>(0)), y(static_cast<T>(0)), z(static_cast<T>(0)), num(0)
{}

// Конструктор по трем координатам и номеру
template<typename T>
template<typename U>
point3_t<T>::point3_t(U n_x, U n_y, U n_z, std::size_t n_num)
    : x(static_cast<T>(n_x)), y(static_cast<T>(n_y)), z(static_cast<T>(n_z)), num(n_num)
{}

// Конструктор из другой точки
template<typename T>
template<typename U>
point3_t<T>::point3_t(const point3_t<U> & p)
{
    x = static_cast<T>(p.x);
    y = static_cast<T>(p.y);
    z = static_cast<T>(p.z);
    num = p.num;
}

// Операторы типа "скобка"
template<typename T>
T & point3_t<T>::operator [] (std::size_t i)
{
    assert(i < 3);
    switch(i)
    {
    case 0:
        return x;
    case 1:
        return y;
    case 2:
        return z;
    }
    return x;
}

template<typename T>
const T & point3_t<T>::operator [] (std::size_t i) const
{
    assert(i < 3);
    switch(i)
    {
    case 0:
        return x;
    case 1:
        return y;
    case 2:
        return z;
    }
    return x;
}

// Оператор меньше (по номеру)
template<typename T>
template<typename U>
bool point3_t<T>::operator < (const point3_t<U> & t) const
{
    return num < t.num;
}

// Оператор равенства (по номеру)
template<typename T>
template<typename U>
bool point3_t<T>::operator == (const point3_t<U> & t) const
{
    return num == t.num;
}

// Оператор присваивания
template<typename T>
point3_t<T> & point3_t<T>::operator = (const point3_t<T> & other)
{
    if(this != & other)
    {
        this->x = other.x;
        this->y = other.y;
        this->z = other.z;
        this->num = other.num;
    }
    return * this;
}

template<typename T>
template<typename U>
point3_t<T> & point3_t<T>::operator = (const point3_t<U> & other)
{
    this->x = static_cast<T>(other.x);
    this->y = static_cast<T>(other.y);
    this->z = static_cast<T>(other.z);
    this->num = other.num;
    return * this;
}

// Вывод
template<typename U>
std::ostream & operator << (std::ostream & os, const point3_t<U> & a)
{
    os << "{ " << a.x << ", " << a.y << ", " << a.z << " }";
    return os;
}

// Проверка, лежит ли точка внутри параллелепипеда (для дерева)
template<typename T>
bool point3_t<T>::inside(double x0, double x1, double y0, double y1, double z0, double z1) const
{
    if(x >= x0 && x <= x1 && y >= y0 && y <= y1 && z >= z0 && z <= z1)
        return true;
    return false;
}

// *************************************************************************************************

}}} // namespace core::containers::geometry

#endif // CONTAINERS_GEOMETRY_POINT3_T_H_INCLUDED
