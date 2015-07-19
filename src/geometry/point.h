#if !defined POINT_H_INCLUDED
#define POINT_H_INCLUDED

#include "../common/common.h"

// Немного typedef'ов
template<typename T> class point_t;
typedef point_t<double> point;
typedef point_t< complex<double> > cpoint;

// Класс точка
template<typename T>
class point_t
{
public:
    // Координаты
    T x, y, z;
    // Номер точки
    size_t num;
    // Конструктор по умолчанию
    point_t();
    // Конструктор по трем координатам
    template<typename U>
    point_t(U x, U y, U z);
    // Конструктор по трем координатам и номеру
    template<typename U>
    point_t(U x, U y, U z, size_t num);
    // Конструктор из другой точки
    template<typename U>
    point_t(const point_t<U> & p);
    // Операторы типа "скобка"
    T & operator [] (size_t i);
    const T & operator [] (size_t i) const;
    // Оператор меньше (по номеру)
    template<typename U>
    bool operator < (const point_t<U> & t) const;
    // Оператор равенства (по номеру)
    template<typename U>
    bool operator == (const point_t<U> & t) const;
    // Оператор присваивания
    point_t<T> & operator = (const point_t<T> & other);
    template<typename U>
    point_t<T> & operator = (const point_t<U> & other);
    // Вывод
    template<typename U>
    friend ostream & operator << (ostream & os, const point_t<U> & a);
    // Проверка, лежит ли точка внутри параллелепипеда (для дерева)
    bool inside(double x0, double x1, double y0, double y1, double z0, double z1) const;
};

// Конструктор по умолчанию
template<typename T>
point_t<T>::point_t()
{
    x = y = z = static_cast<T>(0);
    num = 0;
}

// Конструктор по трем координатам
template<typename T>
template<typename U>
point_t<T>::point_t(U x, U y, U z)
{
    this->x = static_cast<T>(x);
    this->y = static_cast<T>(y);
    this->z = static_cast<T>(z);
    num = 0;
}

// Конструктор по трем координатам и номеру
template<typename T>
template<typename U>
point_t<T>::point_t(U x, U y, U z, size_t num)
{
    this->x = static_cast<T>(x);
    this->y = static_cast<T>(y);
    this->z = static_cast<T>(z);
    this->num = num;
}

// Конструктор из другой точки
template<typename T>
template<typename U>
point_t<T>::point_t(const point_t<U> & p)
{
    x = static_cast<T>(p.x);
    y = static_cast<T>(p.y);
    z = static_cast<T>(p.z);
    num = p.num;
}

// Операторы типа "скобка"
template<typename T>
T & point_t<T>::operator [] (size_t i)
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
const T & point_t<T>::operator [] (size_t i) const
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
bool point_t<T>::operator < (const point_t<U> & t) const
{
    return num < t.num;
}

// Оператор равенства (по номеру)
template<typename T>
template<typename U>
bool point_t<T>::operator == (const point_t<U> & t) const
{
    return num == t.num;
}

// Оператор присваивания
template<typename T>
point_t<T> & point_t<T>::operator = (const point_t<T> & other)
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
point_t<T> & point_t<T>::operator = (const point_t<U> & other)
{
    this->x = static_cast<T>(other.x);
    this->y = static_cast<T>(other.y);
    this->z = static_cast<T>(other.z);
    this->num = other.num;
    return * this;
}

// Вывод
template<typename U>
ostream & operator << (ostream & os, const point_t<U> & a)
{
    os << "{ " << a.x << ", " << a.y << ", " << a.z << " }";
    return os;
}

// Проверка, лежит ли точка внутри параллелепипеда (для дерева)
template<typename T>
bool point_t<T>::inside(double x0, double x1, double y0, double y1, double z0, double z1) const
{
    if(x >= x0 && x <= x1 && y >= y0 && y <= y1 && z >= z0 && z <= z1)
        return true;
    return false;
}

#endif // POINT_H_INCLUDED
