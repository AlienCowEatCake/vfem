#if !defined(CONTAINERS_GEOMETRY_VECTOR3_T_H_INCLUDED)
#define CONTAINERS_GEOMETRY_VECTOR3_T_H_INCLUDED

#include <cstring>
#include <fstream>
#include <complex>

#include "point3_t.h"
#include "../generic/array_t.h"
#include "../generic/matrix_t.h"

namespace core { namespace containers { namespace geometry {

// *************************************************************************************************

// Структура трехмерный вектор
template<typename T>
struct vector3_t
{
    // Элементы
    T x, y, z;
    // Конструктор по-умолчанию
    vector3_t();
    // Конструктор из трех элементов
    template<typename U, typename V, typename R>
    vector3_t(U n_x, V n_y, R n_z);
    // Конструктор из точки
    template<typename U>
    vector3_t(const point3_t<U> & p);
    // Конструктор из двух точек - начала и конца
    template<typename U>
    vector3_t(const point3_t<U> & beg, const point3_t<U> & end);
    // Конструктор из другого вектора
    template<typename U>
    vector3_t(const vector3_t<U> & v);
    // Кастование к точке
    point3_t<T> to_point() const;
    // Операторы типа "скобка"
    T & operator [] (std::size_t i);
    const T & operator [] (std::size_t i) const;
    // Скалярное произведение
    T operator * (const vector3_t<T> & other) const;
    // Векторное произведение
    vector3_t<T> cross(const vector3_t<T> & other) const;
    // Получение сопряженного вектора
    vector3_t<T> conj() const;
    // Норма вектора
    double norm() const;
    // Квадрат нормы вектора
    double norm2() const;
    // Сложение и вычитание векторов
    vector3_t<T> operator + (const vector3_t<T> & other) const;
    vector3_t<T> operator - (const vector3_t<T> & other) const;
    vector3_t<T> & operator += (const vector3_t<T> & other);
    vector3_t<T> & operator -= (const vector3_t<T> & other);
    // Деление вектора на число
    vector3_t<T> operator / (const T & a) const;
    // Умножение вектора на число
    vector3_t<T> operator * (const T & a) const;
    vector3_t<T> & operator *= (const T & a);
    // Умножение числа на вектор
    template<typename U>
    friend vector3_t<U> operator * (const U & a, const vector3_t<U> & vec);
    // Умножение матрицы на вектор
    template<typename U, typename V>
    friend vector3_t<U> operator * (const generic::matrix_t<V, 3, 3> & matr, const vector3_t<U> & vec);
    // Вывод
    template<typename U>
    friend std::ostream & operator << (std::ostream & os, const vector3_t<U> & a);

    // ***** Боль, страдания и унижение :( *****
    // Умножение числа на вектор
    template<typename U, typename V, typename R>
    friend R operator * (const U & a, const vector3_t<V> & vec);
    // Скалярное произведение
    template<typename U, typename V, typename R>
    friend R operator * (const vector3_t<U> & left, const vector3_t<V> & right);
    // Сложение и вычитание векторов
    template<typename U, typename V, typename R>
    friend R operator + (const vector3_t<U> & left, const vector3_t<V> & right);
    template<typename U, typename V, typename R>
    friend R operator - (const vector3_t<U> & left, const vector3_t<V> & right);
};

// *************************************************************************************************

// Конструктор по-умолчанию
template<typename T>
vector3_t<T>::vector3_t()
    : x(static_cast<T>(0)), y(static_cast<T>(0)), z(static_cast<T>(0))
{}

// Конструктор из трех элементов
template<typename T>
template<typename U, typename V, typename R>
vector3_t<T>::vector3_t(U n_x, V n_y, R n_z)
    : x(static_cast<T>(n_x)), y(static_cast<T>(n_y)), z(static_cast<T>(n_z))
{}

// Конструктор из точки
template<typename T>
template<typename U>
vector3_t<T>::vector3_t(const point3_t<U> & p)
    : x(static_cast<T>(p.x)), y(static_cast<T>(p.y)), z(static_cast<T>(p.z))
{}

// Конструктор из двух точек - начала и конца
template<typename T>
template<typename U>
vector3_t<T>::vector3_t(const point3_t<U> & beg, const point3_t<U> & end)
    : x(static_cast<T>(end.x - beg.x)), y(static_cast<T>(end.y - beg.y)), z(static_cast<T>(end.z - beg.z))
{}

// Конструктор из другого вектора
template<typename T>
template<typename U>
vector3_t<T>::vector3_t(const vector3_t<U> & v)
    : x(static_cast<T>(v.x)), y(static_cast<T>(v.y)), z(static_cast<T>(v.z))
{}

// Кастование к точке
template<typename T>
point3_t<T> vector3_t<T>::to_point() const
{
    return point3_t<T>(x, y, z);
}

// Операторы типа "скобка"
template<typename T>
T & vector3_t<T>::operator [] (std::size_t i)
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
const T & vector3_t<T>::operator [] (std::size_t i) const
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

// Скалярное произведение
template<typename T>
T vector3_t<T>::operator * (const vector3_t<T> & other) const
{
    return x * other.x + y * other.y + z * other.z;
}

// Векторное произведение
template<typename T>
vector3_t<T> vector3_t<T>::cross(const vector3_t<T> & other) const
{
    return vector3_t<T>(y * other.z - z * other.y,
                        z * other.x - x * other.z,
                        x * other.y - y * other.x);
}

namespace {

// Получение сопряженного вектора (комплексного)
template<typename T>
vector3_t< std::complex<T> > vector3_t_conj(const vector3_t< std::complex<T> > & vec)
{
    return vector3_t< std::complex<T> >(std::conj(vec.x), std::conj(vec.y), std::conj(vec.z));
}

// Получение сопряженного вектора (действительного)
template<typename T>
vector3_t<T> vector3_t_conj(const vector3_t<T> & vec)
{
    return vec;
}

} // namespace

// Получение сопряженного вектора
template<typename T>
vector3_t<T> vector3_t<T>::conj() const
{
    return vector3_t_conj(* this);
}
// Норма вектора
template<typename T>
double vector3_t<T>::norm() const
{
    return sqrt(norm2());
}

namespace {

// Квадрат нормы вектора (комплексного)
template<typename T>
T vector3_t_norm2(const vector3_t<std::complex<T> > & v)
{
    T x_re = v.x.real(), x_im = v.x.imag();
    T y_re = v.y.real(), y_im = v.y.imag();
    T z_re = v.z.real(), z_im = v.z.imag();
    return x_re * x_re + x_im * x_im +
           y_re * y_re + y_im * y_im +
           z_re * z_re + z_im * z_im;
}

// Квадрат нормы вектора (действительного)
template<typename T>
T vector3_t_norm2(const vector3_t<T> & v)
{
    return v.x * v.x + v.y * v.y + v.z * v.z;
}

} // namespace

// Квадрат нормы вектора
template<typename T>
double vector3_t<T>::norm2() const
{
    return vector3_t_norm2(* this);
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

template<typename T>
vector3_t<T> & vector3_t<T>::operator += (const vector3_t<T> & other)
{
    x += other.x;
    y += other.y;
    z += other.z;
    return * this;
}

template<typename T>
vector3_t<T> & vector3_t<T>::operator -= (const vector3_t<T> & other)
{
    x -= other.x;
    y -= other.y;
    z -= other.z;
    return * this;
}

// Деление вектора на число
template<typename T>
vector3_t<T> vector3_t<T>::operator / (const T & a) const
{
    return vector3_t<T>(x / a, y / a, z / a);
}

// Умножение вектора на число
template<typename T>
vector3_t<T> vector3_t<T>::operator * (const T & a) const
{
    return vector3_t<T>(x * a, y * a, z * a);
}


template<typename T>
vector3_t<T> & vector3_t<T>::operator *= (const T & a)
{
    x *= a;
    y *= a;
    z *= a;
    return * this;
}

// Умножение числа на вектор
template<typename U>
vector3_t<U> operator * (const U & a, const vector3_t<U> & vec)
{
    return vector3_t<U>(a * vec.x, a * vec.y, a * vec.z);
}

// Умножение матрицы на вектор
template<typename U, typename V>
vector3_t<U> operator * (const generic::matrix_t<V, 3, 3> & matr, const vector3_t<U> & vec)
{
    vector3_t<U> result(static_cast<U>(0), static_cast<U>(0), static_cast<U>(0));
    for(std::size_t i = 0; i < 3; i++)
        for(std::size_t j = 0; j < 3; j++)
            result[i] += static_cast<U>(matr[i][j]) * vec[j];
    return result;
}

// Вывод
template<typename U>
std::ostream & operator << (std::ostream & os, const vector3_t<U> & a)
{
    os << "{ " << a.x << ", " << a.y << ", " << a.z << " }";
    return os;
}

// *************************************************************************************************

// ***** Боль, страдания и унижение :( *****

// Умножение числа на вектор
inline vector3_t<double> operator * (const float & a, const vector3_t<double> & vec)
{
    return vector3_t<double>(vec.x * static_cast<double>(a),
                             vec.y * static_cast<double>(a),
                             vec.z * static_cast<double>(a));
}

inline vector3_t< std::complex<float> > operator * (const float & a, const vector3_t< std::complex<float> > & vec)
{
    return vector3_t< std::complex<float> >(vec.x * a,
                                            vec.y * a,
                                            vec.z * a);
}

inline vector3_t< std::complex<double> > operator * (const float & a, const vector3_t< std::complex<double> > & vec)
{
    return vector3_t< std::complex<double> >(vec.x * static_cast<double>(a),
                                             vec.y * static_cast<double>(a),
                                             vec.z * static_cast<double>(a));
}

inline vector3_t<float> operator * (const double & a, const vector3_t<float> & vec)
{
    return vector3_t<float>(vec.x * static_cast<float>(a),
                            vec.y * static_cast<float>(a),
                            vec.z * static_cast<float>(a));
}

inline vector3_t< std::complex<float> > operator * (const double & a, const vector3_t< std::complex<float> > & vec)
{
    return vector3_t< std::complex<float> >(vec.x * static_cast<float>(a),
                                            vec.y * static_cast<float>(a),
                                            vec.z * static_cast<float>(a));
}

inline vector3_t< std::complex<double> > operator * (const double & a, const vector3_t< std::complex<double> > & vec)
{
    return vector3_t< std::complex<double> >(vec.x * a,
                                             vec.y * a,
                                             vec.z * a);
}

inline vector3_t< std::complex<float> > operator * (const std::complex<float> & a, const vector3_t<float> & vec)
{
    return vector3_t< std::complex<float> >(vec.x * a,
                                            vec.y * a,
                                            vec.z * a);
}

inline vector3_t< std::complex<double> > operator * (const std::complex<float> & a, const vector3_t<double> & vec)
{
    return vector3_t< std::complex<double> >(vec.x * std::complex<double>(static_cast<double>(a.real()), static_cast<double>(a.imag())),
                                             vec.y * std::complex<double>(static_cast<double>(a.real()), static_cast<double>(a.imag())),
                                             vec.z * std::complex<double>(static_cast<double>(a.real()), static_cast<double>(a.imag())));
}

inline vector3_t< std::complex<double> > operator * (const std::complex<float> & a, const vector3_t< std::complex<double> > & vec)
{
    return vector3_t< std::complex<double> >(vec.x * std::complex<double>(static_cast<double>(a.real()), static_cast<double>(a.imag())),
                                             vec.y * std::complex<double>(static_cast<double>(a.real()), static_cast<double>(a.imag())),
                                             vec.z * std::complex<double>(static_cast<double>(a.real()), static_cast<double>(a.imag())));
}

inline vector3_t< std::complex<float> > operator * (const std::complex<double> & a, const vector3_t<float> & vec)
{
    return vector3_t< std::complex<float> >(vec.x * std::complex<float>(static_cast<float>(a.real()), static_cast<float>(a.imag())),
                                            vec.y * std::complex<float>(static_cast<float>(a.real()), static_cast<float>(a.imag())),
                                            vec.z * std::complex<float>(static_cast<float>(a.real()), static_cast<float>(a.imag())));
}

inline vector3_t< std::complex<double> > operator * (const std::complex<double> & a, const vector3_t<double> & vec)
{
    return vector3_t< std::complex<double> >(vec.x * a,
                                             vec.y * a,
                                             vec.z * a);
}

inline vector3_t< std::complex<float> > operator * (const std::complex<double> & a, const vector3_t< std::complex<float> > & vec)
{
    return vector3_t< std::complex<float> >(vec.x * std::complex<float>(static_cast<float>(a.real()), static_cast<float>(a.imag())),
                                            vec.y * std::complex<float>(static_cast<float>(a.real()), static_cast<float>(a.imag())),
                                            vec.z * std::complex<float>(static_cast<float>(a.real()), static_cast<float>(a.imag())));
}

// Скалярное произведение
inline double operator * (const vector3_t<float> & left, const vector3_t<double> & right)
{
    return vector3_t<double>(left) * right;
}

inline std::complex<float> operator * (const vector3_t<float> & left, const vector3_t< std::complex<float> > & right)
{
    return vector3_t< std::complex<float> >(left) * right;
}

inline std::complex<double> operator * (const vector3_t<float> & left, const vector3_t< std::complex<double> > & right)
{
    return vector3_t< std::complex<double> >(vector3_t<double>(left)) * right;
}

inline double operator * (const vector3_t<double> & left, const vector3_t<float> & right)
{
    return left * vector3_t<double>(right);
}

inline std::complex<float> operator * (const vector3_t<double> & left, const vector3_t< std::complex<float> > & right)
{
    return vector3_t< std::complex<float> >(vector3_t<float>(left)) * right;
}

inline std::complex<double> operator * (const vector3_t<double> & left, const vector3_t< std::complex<double> > & right)
{
    return vector3_t< std::complex<double> >(left) * right;
}

inline std::complex<float> operator * (const vector3_t< std::complex<float> > & left, const vector3_t<float> & right)
{
    return left * vector3_t< std::complex<float> >(right);
}

inline std::complex<float> operator * (const vector3_t< std::complex<float> > & left, const vector3_t<double> & right)
{
    return left * vector3_t< std::complex<float> >(vector3_t<float>(right));
}

inline std::complex<double> operator * (const vector3_t< std::complex<float> > & left, const vector3_t< std::complex<double> > & right)
{
    return vector3_t< std::complex<double> >(left) * right;
}

inline std::complex<double> operator * (const vector3_t< std::complex<double> > & left, const vector3_t<float> & right)
{
    return left * vector3_t< std::complex<double> >(vector3_t<double>(right));
}

inline std::complex<double> operator * (const vector3_t< std::complex<double> > & left, const vector3_t<double> & right)
{
    return left * vector3_t< std::complex<double> >(right);
}

inline std::complex<double> operator * (const vector3_t< std::complex<double> > & left, const vector3_t< std::complex<float> > & right)
{
    return left * vector3_t< std::complex<double> >(right);
}

// Сложение векторов
inline vector3_t<double> operator + (const vector3_t<float> & left, const vector3_t<double> & right)
{
    return vector3_t<double>(left) + right;
}

inline vector3_t< std::complex<float> > operator + (const vector3_t<float> & left, const vector3_t< std::complex<float> > & right)
{
    return vector3_t< std::complex<float> >(left) + right;
}

inline vector3_t< std::complex<double> > operator + (const vector3_t<float> & left, const vector3_t< std::complex<double> > & right)
{
    return vector3_t< std::complex<double> >(vector3_t<double>(left)) + right;
}

inline vector3_t<double> operator + (const vector3_t<double> & left, const vector3_t<float> & right)
{
    return left + vector3_t<double>(right);
}

inline vector3_t< std::complex<float> > operator + (const vector3_t<double> & left, const vector3_t< std::complex<float> > & right)
{
    return vector3_t< std::complex<float> >(vector3_t<float>(left)) + right;
}

inline vector3_t< std::complex<double> > operator + (const vector3_t<double> & left, const vector3_t< std::complex<double> > & right)
{
    return vector3_t< std::complex<double> >(left) + right;
}

inline vector3_t< std::complex<float> > operator + (const vector3_t< std::complex<float> > & left, const vector3_t<float> & right)
{
    return left + vector3_t< std::complex<float> >(right);
}

inline vector3_t< std::complex<float> > operator + (const vector3_t< std::complex<float> > & left, const vector3_t<double> & right)
{
    return left + vector3_t< std::complex<float> >(vector3_t<float>(right));
}

inline vector3_t< std::complex<double> > operator + (const vector3_t< std::complex<float> > & left, const vector3_t< std::complex<double> > & right)
{
    return vector3_t< std::complex<double> >(left) + right;
}

inline vector3_t< std::complex<double> > operator + (const vector3_t< std::complex<double> > & left, const vector3_t<float> & right)
{
    return left + vector3_t< std::complex<double> >(vector3_t<double>(right));
}

inline vector3_t< std::complex<double> > operator + (const vector3_t< std::complex<double> > & left, const vector3_t<double> & right)
{
    return left + vector3_t< std::complex<double> >(right);
}

inline vector3_t< std::complex<double> > operator + (const vector3_t< std::complex<double> > & left, const vector3_t< std::complex<float> > & right)
{
    return left + vector3_t< std::complex<double> >(right);
}

// Вычитание векторов
inline vector3_t<double> operator - (const vector3_t<float> & left, const vector3_t<double> & right)
{
    return vector3_t<double>(left) - right;
}

inline vector3_t< std::complex<float> > operator - (const vector3_t<float> & left, const vector3_t< std::complex<float> > & right)
{
    return vector3_t< std::complex<float> >(left) - right;
}

inline vector3_t< std::complex<double> > operator - (const vector3_t<float> & left, const vector3_t< std::complex<double> > & right)
{
    return vector3_t< std::complex<double> >(vector3_t<double>(left)) - right;
}

inline vector3_t<double> operator - (const vector3_t<double> & left, const vector3_t<float> & right)
{
    return left - vector3_t<double>(right);
}

inline vector3_t< std::complex<float> > operator - (const vector3_t<double> & left, const vector3_t< std::complex<float> > & right)
{
    return vector3_t< std::complex<float> >(vector3_t<float>(left)) - right;
}

inline vector3_t< std::complex<double> > operator - (const vector3_t<double> & left, const vector3_t< std::complex<double> > & right)
{
    return vector3_t< std::complex<double> >(left) - right;
}

inline vector3_t< std::complex<float> > operator - (const vector3_t< std::complex<float> > & left, const vector3_t<float> & right)
{
    return left - vector3_t< std::complex<float> >(right);
}

inline vector3_t< std::complex<float> > operator - (const vector3_t< std::complex<float> > & left, const vector3_t<double> & right)
{
    return left - vector3_t< std::complex<float> >(vector3_t<float>(right));
}

inline vector3_t< std::complex<double> > operator - (const vector3_t< std::complex<float> > & left, const vector3_t< std::complex<double> > & right)
{
    return vector3_t< std::complex<double> >(left) - right;
}

inline vector3_t< std::complex<double> > operator - (const vector3_t< std::complex<double> > & left, const vector3_t<float> & right)
{
    return left - vector3_t< std::complex<double> >(vector3_t<double>(right));
}

inline vector3_t< std::complex<double> > operator - (const vector3_t< std::complex<double> > & left, const vector3_t<double> & right)
{
    return left - vector3_t< std::complex<double> >(right);
}

inline vector3_t< std::complex<double> > operator - (const vector3_t< std::complex<double> > & left, const vector3_t< std::complex<float> > & right)
{
    return left - vector3_t< std::complex<double> >(right);
}

// *************************************************************************************************

}}} // namespace core::containers::geometry

#endif // CONTAINERS_GEOMETRY_VECTOR3_T_H_INCLUDED
