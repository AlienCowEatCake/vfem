#if !defined VECTOR3_H_INCLUDED
#define VECTOR3_H_INCLUDED

#include "../common/common.h"
#include "../common/matrix.h"
#include "../geometry/point.h"

// Немного typedef'ов
template<typename T> class vector3_t;
typedef vector3_t<double> vector3;
typedef vector3_t< complex<double> > cvector3;

// Класс трехмерный вектор
template<typename T>
class vector3_t
{
public:
    // Элементы
    T x, y, z;
    // Констуктор по умолчанию
    vector3_t();
    // Конструктор из трех элементов
    template<typename U, typename V, typename R>
    vector3_t(U x, V y, R z);
    // Конструктор из точки
    template<typename U>
    vector3_t(const point_t<U> & p);
    // Конструктор из двух точек - начала и конца
    template<typename U>
    vector3_t(const point_t<U> & beg_, const point_t<U> & end_);
    // Конструктор из другого вектора
    template<typename U>
    vector3_t(const vector3_t<U> & v);
    // Кастование к точке
    point_t<T> pnt() const;
    // Норма вектора
    double norm() const;
    // Векторное произведение
    vector3_t<T> cross(const vector3_t<T> & other) const;
    // Операторы типа "скобка"
    T & operator [] (size_t i);
    T operator [] (size_t i) const;
    // Скалярное произведение
    T operator * (const vector3_t<T> & other) const;
    template<typename U, typename V>
    V operator * (const vector3_t<U> & other) const;
    // Сложение и вычитание векторов
    vector3_t<T> operator + (const vector3_t<T> & other) const;
    vector3_t<T> operator - (const vector3_t<T> & other) const;
    // Деление вектора на число
    vector3_t<T> operator / (const T & a) const;
    // Умножение числа на вектор
    template<typename U>
    friend vector3_t<U> operator * (const U & a, const vector3_t<U> & vec);
    template<typename U, typename V, typename R>
    friend vector3_t<R> operator * (const U & a, const vector3_t<V> & vec);
    // Умножение матрицы на вектор
    template<typename U, typename V>
    friend vector3_t<U> operator * (const matrix_t<V, 3, 3> & matr, const vector3_t<U> & vec);
    // Вывод
    template<typename U>
    friend ostream & operator << (ostream & os, const vector3_t<U> & a);
    // Получение сопряженного вектора
    vector3_t<T> cj() const;
};

// Констуктор по умолчанию
template<typename T>
vector3_t<T>::vector3_t()
{
    x = y = z = static_cast<T>(0);
}

// Конструктор из трех элементов
template<typename T>
template<typename U, typename V, typename R>
vector3_t<T>::vector3_t(U x, V y, R z)
{
    this->x = static_cast<T>(x);
    this->y = static_cast<T>(y);
    this->z = static_cast<T>(z);
}

// Конструктор из точки
template<typename T>
template<typename U>
vector3_t<T>::vector3_t(const point_t<U> & p)
{
    x = static_cast<T>(p.x);
    y = static_cast<T>(p.y);
    z = static_cast<T>(p.z);
}

// Конструктор из двух точек - начала и конца
template<typename T>
template<typename U>
vector3_t<T>::vector3_t(const point_t<U> & beg_, const point_t<U> & end_)
{
    x = static_cast<T>(end_.x - beg_.x);
    y = static_cast<T>(end_.y - beg_.y);
    z = static_cast<T>(end_.z - beg_.z);
}

// Конструктор из другого вектора
template<typename T>
template<typename U>
vector3_t<T>::vector3_t(const vector3_t<U> & v)
{
    x = static_cast<T>(v.x);
    y = static_cast<T>(v.y);
    z = static_cast<T>(v.z);
}

// Кастование к точке
template<typename T>
point_t<T> vector3_t<T>::pnt() const
{
    return point_t<T>(x, y, z);
}

// Векторное произведение
template<typename T>
vector3_t<T> vector3_t<T>::cross(const vector3_t<T> & other) const
{
    return vector3_t<T>(y * other.z - z * other.y,
                      z * other.x - x * other.z,
                      x * other.y - y * other.x);
}

// Операторы типа "скобка"
template<typename T>
T & vector3_t<T>::operator [] (size_t i)
{
    switch(i)
    {
    case 0:
        return x;
    case 1:
        return y;
    case 2:
        return z;
    default:
        cerr << "Error: Unknown index " << i << " at vector " << * this << endl;
        throw ADDRESSING_ERROR;
    }
}

template<typename T>
T vector3_t<T>::operator [] (size_t i) const
{
    switch(i)
    {
    case 0:
        return x;
    case 1:
        return y;
    case 2:
        return z;
    default:
        cerr << "Error: Unknown index " << i << " at vector " << * this << endl;
        throw ADDRESSING_ERROR;
    }
}

// Скалярное произведение
template<typename T>
T vector3_t<T>::operator * (const vector3_t<T> & other) const
{
    return x * other.x + y * other.y + z * other.z;
}

// Сложение и вычитание векторов
template<typename T>
vector3_t<T> vector3_t<T>::operator + (const vector3_t<T> & other) const
{
    return vector3_t<T>(x + other.x, y + other.y, z + other.z);
}

template<typename T>
vector3_t<T> vector3_t<T>::operator - (const vector3_t<T> & other) const
{
    return vector3_t<T>(x - other.x, y - other.y, z - other.z);
}

// Деление вектора на число
template<typename T>
vector3_t<T> vector3_t<T>::operator / (const T & a) const
{
    return vector3_t<T>(x / a, y / a, z / a);
}

// Умножение числа на вектор
template<typename U>
vector3_t<U> operator * (const U & a, const vector3_t<U> & vec)
{
    return vector3_t<U>(a * vec.x, a * vec.y, a * vec.z);
}

// Умножение матрицы на вектор
template<typename U, typename V>
vector3_t<U> operator * (const matrix_t<V, 3, 3> & matr, const vector3_t<U> & vec)
{
    vector3_t<U> result(static_cast<U>(0), static_cast<U>(0), static_cast<U>(0));
    for(size_t i = 0; i < 3; i++)
        for(size_t j = 0; j < 3; j++)
            result[i] += static_cast<U>(matr[i][j]) * vec[j];
    return result;
}

// Вывод
template<typename U>
ostream & operator << (ostream & os, const vector3_t<U> & a)
{
    os << "{ " << a.x << ", " << a.y << ", " << a.z << " }";
    return os;
}

// Норма вектора (действительная)
template<>
double vector3_t<double>::norm() const;

// Норма вектора (комплексная)
template<>
double vector3_t< complex<double> >::norm() const;

// Скалярное произведение
template<> template<>
complex<double> vector3_t<double>::operator * (const vector3_t< complex<double> > & other) const;

template<> template<>
complex<double> vector3_t< complex<double> >::operator * (const vector3_t<double> & other) const;

// Умножение числа на вектор
vector3_t< complex<double> > operator * (const double & a, const vector3_t< complex<double> > & vec);

vector3_t< complex<double> > operator * (const complex<double> & a, const vector3_t<double> & vec);

// Получение сопряженного вектора
template<>
vector3_t< complex<double> > vector3_t< complex<double> >::cj() const;

#endif // VECTOR3_H_INCLUDED
