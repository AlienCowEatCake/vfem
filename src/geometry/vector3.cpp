#include "vector3.h"

// Норма вектора (действительная)
template<>
double vector3_t<double>::norm() const
{
    return sqrt(x * x + y * y + z * z);
}

// Норма вектора (комплексная)
template<>
double vector3_t< complex<double> >::norm() const
{
    return sqrt((conj(x) * x + conj(y) * y + conj(z) * z).real());
}

// Скалярное произведение
template<> template<>
complex<double> vector3_t<double>::operator * (const vector3_t< complex<double> > & other) const
{
    return complex<double>(x, 0.0) * other.x +
           complex<double>(y, 0.0) * other.y +
           complex<double>(z, 0.0) * other.z;
}

template<> template<>
complex<double> vector3_t< complex<double> >::operator * (const vector3_t<double> & other) const
{
    return x * complex<double>(other.x, 0.0) +
           y * complex<double>(other.y, 0.0) +
           z * complex<double>(other.z, 0.0);
}

// Умножение числа на вектор
vector3_t< complex<double> > operator * (const double & a, const vector3_t< complex<double> > & vec)
{
    return vector3_t< complex<double> >(a * vec.x, a * vec.y, a * vec.z);
}

vector3_t< complex<double> > operator * (const complex<double> & a, const vector3_t<double> & vec)
{
    return vector3_t< complex<double> >(a * vec.x, a * vec.y, a * vec.z);
}

// Получение сопряженного вектора
template<>
vector3_t< complex<double> > vector3_t< complex<double> >::cj() const
{
    return vector3_t< complex<double> >(conj(x), conj(y), conj(z));
}
